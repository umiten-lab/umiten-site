#!/usr/bin/env python3
"""海点: 沖の波(Open-Meteoの波浪モデル)と風(気象庁モデル)から、各ポイントの岸での波のサイズを計算し、
docs/data/forecast.json に書き出す。標準ライブラリだけで動く。"""
import json, math, os, ssl, sys, time, datetime, pathlib, urllib.parse, urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
G = 9.81
JST = datetime.timezone(datetime.timedelta(hours=9))
LAG_UP = 3.0     # 面が悪→良に向かう時定数(h)。荒れた水面は、風が変わってもすぐには落ち着かない
LAG_DOWN = 2.0   # 面が良→悪に向かう時定数(h)。ずっとオフショアなら、風が変わってもすぐには崩れない
HOLD_BAD, HOLD_GOOD = 2, 3   # 風が悪化(改善)したと数えるまでに、その状態が続くべき時間(h)
CTX = ssl._create_unverified_context() if os.environ.get('UMITEN_INSECURE') else None   # 証明書の検証を切る(試用のみ)
MARINE = ("wave_height,wave_period,wave_direction,wind_wave_height,wind_wave_period,wind_wave_direction,"
          "swell_wave_height,swell_wave_period,swell_wave_direction,secondary_swell_wave_height,"
          "secondary_swell_wave_period,secondary_swell_wave_direction,sea_surface_temperature,sea_level_height_msl")

def get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "umiten-forecast/0.1"})
            with urllib.request.urlopen(req, timeout=30, context=CTX) as r:
                return json.load(r)
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))

def angdiff(a, b):
    return abs((a - b + 180) % 360 - 180)

def wavelength(T, h):
    L0 = G * T * T / (2 * math.pi); L = L0
    for _ in range(60):
        L = L0 * math.tanh(2 * math.pi * h / L)
    return L

def to_shore(H, T, theta, facing, depth):
    """沖の波(高さH, 周期T, 来る向きtheta)を、水深depthの岸近くの波高に直す。
    線形波理論: スネルの法則による屈折(Kr)と浅水変形(Ks)。海底の等深線は岸と平行と仮定。"""
    if not H or not T or T < 1 or theta is None:
        return 0.0
    a0 = angdiff(theta, facing)
    if a0 >= 89:
        return 0.0
    L0 = G * T * T / (2 * math.pi); L = wavelength(T, depth)
    s = (L / L0) * math.sin(math.radians(a0))
    if s >= 1:
        return 0.0
    kr = math.sqrt(math.cos(math.radians(a0)) / math.sqrt(1 - s * s))
    k = 2 * math.pi / L
    n = 0.5 * (1 + 2 * k * depth / math.sinh(2 * k * depth))
    ks = math.sqrt(0.5 * L0 / (n * L))
    return H * kr * ks

def wind_kind(frm, facing):
    """風の吹いてくる向きと海の向きの差で、Surf life と同じ5段階に分ける。"""
    if frm is None:
        return None
    d = angdiff(frm, facing)
    return "on" if d <= 30 else "sideon" if d <= 60 else "side" if d < 120 else "sideoff" if d < 150 else "off"

def wind_quality(spd, frm, facing):
    """風だけから見た面の良さ(0〜1)。オフショアで弱いほど高く、オンショアほど低い。"""
    kind = wind_kind(frm, facing)
    if spd is None or kind is None:
        return None
    t = {"off": [(6, 1.0), (8, .85), (10, .6), (13, .35)], "sideoff": [(5, .9), (7, .75), (10, .5), (13, .25)],
         "side": [(3, .65), (5, .45), (8, .25)], "sideon": [(2, .35), (4, .18)], "on": [(1.5, .2), (3, .08)]}[kind]
    q = next((v for lim, v in t if spd <= lim), .1 if kind in ("off", "sideoff", "side") else 0.0)
    if spd <= 1.0:
        return 0.85  # ほぼ無風: 風向きによる優劣をつけない
    if spd <= 1.5:
        return (1.5 - spd) * 2 * 0.85 + (spd - 1.0) * 2 * max(q, 0.6)
    return q

def surface_quality(wqs, up=LAG_UP, down=LAG_DOWN):
    """風の良さ(wqs: 1時間ごと)から、水面の良さを時間差つきで求める。
    風が変わっても、悪化は HOLD_BAD 時間、改善は HOLD_GOOD 時間続いてから反映を始め、その後も指数的にゆっくり変わる。"""
    f, last = [], None
    for v in wqs:                                   # 欠測は直前の値で埋める
        last = v if v is not None else last
        f.append(last)
    out, prev = [], None
    for j, cur in enumerate(f):
        if cur is None:
            out.append(None); continue
        prev = cur if prev is None else prev
        bad = max(f[max(0, j - k)] for k in range(HOLD_BAD))     # 直近 HOLD_BAD 時間の最良
        good = min(f[max(0, j - k)] for k in range(HOLD_GOOD))   # 直近 HOLD_GOOD 時間の最悪
        if bad < prev:
            prev += (bad - prev) * (1 - math.exp(-1 / down))
        elif good > prev:
            prev += (good - prev) * (1 - math.exp(-1 / up))
        out.append(prev)
    return out

def wind_score_cap(kind, spd):
    """現在の風向・風速から、その時間帯のスコア上限を設定する。

    面の良さ Q は過去の風を考慮するが、現在オンショアに変わった場合に
    古い良好な Q だけで高得点にならないようにする。試験的な経験則。
    """
    if kind is None or spd is None:
        return None
    if spd <= 1.5:  # ほぼ無風。風向きの影響は小さい
        return 100
    base = {"off": 100, "sideoff": 82, "side": 58, "sideon": 42, "on": 35}[kind]
    if spd < 3.0:  # 1.5〜3 m/s の間は方位による上限を徐々に効かせる
        return round(100 - (100 - base) * (spd - 1.5) / 1.5)
    # 風速が強くなるほど減点。オンショア・サイドオンを特に厳しくする。
    start, loss, floor = {
        "off": (8.0, 4.0, 65),
        "sideoff": (6.0, 4.0, 50),
        "side": (4.0, 4.0, 30),
        "sideon": (3.0, 5.0, 15),
        "on": (2.0, 6.0, 10),
    }[kind]
    return max(floor, round(base - loss * max(0.0, spd - start)))


# 日本式の波サイズ表記(当サイトの簡易岸波高に対応させた暫定値)
# ショートは胸〜肩(中心1.35 m)、ロングは腹〜胸(中心1.08 m)。
BOARD_PROFILES = {
    "short": {"mu": 1.35, "sigma": 0.48},
    "long":  {"mu": 1.08, "sigma": 0.38},
}


# 風向の正規分布型評価（角度はオフショア中心からのずれ）。
# ロングは鋭いピーク、ショートは広めのピークを与える。
# 数値は海点の暫定的なサーフ適性パラメータで、観測で校正したものではない。
WIND_DIRECTION_SIGMA = {"short": 80.0, "long": 60.0}  # degree; ショートは広め、ロングはやや鋭め


def offshore_suitability(frm, facing, spd, board="short"):
    """オフショア風向を中心としたガウス型の風向適性 (0–1)。

    方位は風の「吹いてくる向き」。岸が海に向いている facing に対し、
    最良のオフショアは facing+180 度。弱風（≦1m/s）では
    風向の不確実性/実際の海面影響が大きいので角度減点しない。
    3m/s 以上では板種別のガウス型ペナルティを十分に効かせる。
    """
    if frm is None or facing is None or spd is None:
        return None
    theta = angdiff(frm, (facing + 180) % 360)
    sigma = WIND_DIRECTION_SIGMA[board]
    gaussian = math.exp(-0.5 * (theta / sigma) ** 2)
    # 完全なゼロにせず、微弱/瞬間的な方位変化で評価が壊れないようにする。
    floor = 0.18 if board == "long" else 0.20
    strength = max(0.0, min(1.0, (spd - 1.0) / 2.0))
    return 1.0 - strength * (1.0 - (floor + (1.0 - floor) * gaussian))


def size_suitability(a, board="short"):
    """正規分布の指数部を利用したサイズ適合度 (0–1)。統計的な確率ではない。"""
    if a is None:
        return None
    profile = BOARD_PROFILES[board]
    return math.exp(-0.5 * ((max(0.0, a) - profile["mu"]) / profile["sigma"]) ** 2)


def score(a, per, q, kind=None, spd=None, board="short", wind_dir=None, facing=None):
    """板種別の波予報スコア（暫定）：面/風70・サイズ20・周期10。

    サイズ適合度（m）とオフショア中心の風向適合度（度）に
    別々のガウス関数を用いる。ロングは風向分布を鋭く、
    ショートはやや広くし、ロングの横風評価を過度に厳しくしない。
    """
    if a is None or q is None or kind is None or spd is None:
        return None
    a = max(0.0, a)
    q = max(0.0, min(1.0, q))
    fit = size_suitability(a, board)
    period_points = 10 if per >= 10 else 7 if per >= 8 else 3 if per >= 6 else 0
    base = 70 * q + 20 * fit + period_points
    raw = base * (0.35 + 0.65 * fit)

    # 加算式での過大評価を避け、風向の「角度」も連続的に評価する。
    # ほぼ無風では方位を優劣の根拠にしない。
    alignment = offshore_suitability(wind_dir, facing, spd, board)
    if alignment is not None:
        raw *= alignment

    size_cap = (25 if a < 0.35 else
                35 if a < 0.50 else
                39 if a < 0.70 else
                69 if a < 0.90 else
                83 if a < 1.00 else 100)
    period_cap = 59 if per < 6 else 79 if per < 8 else 100
    wind_cap = wind_score_cap(kind, spd)
    return max(0, min(100, round(min(raw, size_cap, period_cap, wind_cap))))


def source_grid(payload):
    """Open-Meteo API の応答に含まれる代表格子点座標。観測所の位置ではない。"""
    lat, lon = payload.get("latitude"), payload.get("longitude")
    return {"lat": round(float(lat), 5) if lat is not None else None,
            "lon": round(float(lon), 5) if lon is not None else None}


def forecast(p):
    q = urllib.parse.urlencode({"latitude": p["lat"], "longitude": p["lon"], "timezone": "Asia/Tokyo", "forecast_days": 7})
    marine_response = get(f"https://marine-api.open-meteo.com/v1/marine?{q}&hourly={MARINE}")
    m = marine_response["hourly"]
    wq = f"https://api.open-meteo.com/v1/forecast?{q}&hourly=wind_speed_10m,wind_direction_10m&wind_speed_unit=ms&past_days=1"   # 前日分は、風向きの変化の履歴(面のラグ)の計算用
    try:
        wind_response = get(wq + "&models=jma_seamless")
        w = wind_response["hourly"]
        assert any(v is not None for v in w["wind_speed_10m"])
        wind_model = "Open-Meteo / JMA Seamless（気象庁モデル）"
    except Exception:
        wind_response = get(wq)
        w = wind_response["hourly"]
        wind_model = "Open-Meteo / 自動選択モデル（JMA優先取得失敗時）"
    wi = {t: i for i, t in enumerate(w["time"])}
    try:   # 天気と気温(表示用。取れなくても続行)
        weather_response = get(f"https://api.open-meteo.com/v1/forecast?{q}&hourly=temperature_2m,weather_code")
        x = weather_response["hourly"]
        xi = {t: i for i, t in enumerate(x["time"])}
    except Exception:
        weather_response, x, xi = None, None, {}
    facing, depth, scale = p["facing"], p.get("depth_m", 4.0), p.get("scale", 1.0)
    up, down = p.get("lag_up_h", LAG_UP), p.get("lag_down_h", LAG_DOWN)
    wqs = [wind_quality(w["wind_speed_10m"][j], w["wind_direction_10m"][j], facing) for j in range(len(w["time"]))]
    Q = dict(zip(w["time"], surface_quality(wqs, up, down)))   # 風の履歴から見た面の良さ
    today = datetime.datetime.now(JST).strftime("%Y-%m-%d")
    rows = []
    for i, t in enumerate(m["time"]):
        if t[:10] < today:
            continue
        g = lambda k: m[k][i]
        comps = [(g("swell_wave_height"), g("swell_wave_period"), g("swell_wave_direction")),
                 (g("secondary_swell_wave_height"), g("secondary_swell_wave_period"), g("secondary_swell_wave_direction")),
                 (g("wind_wave_height"), g("wind_wave_period"), g("wind_wave_direction"))]
        comps = [c for c in comps if c[0]]
        tot = g("wave_height") or 0
        rest = math.sqrt(max(0.0, tot * tot - sum(c[0] ** 2 for c in comps)))   # 上の3成分に入らない残り
        if rest > 0.05:
            comps.append((rest, g("wave_period"), g("wave_direction")))
        shore = [(to_shore(h, T, d, facing, depth), T, d) for h, T, d in comps]
        a = scale * math.sqrt(sum(s[0] ** 2 for s in shore)) if shore else 0.0
        top = max(shore, key=lambda s: s[0]) if shore else (0, g("wave_period"), g("wave_direction"))
        j = wi.get(t)
        spd = w["wind_speed_10m"][j] if j is not None else None
        wd = w["wind_direction_10m"][j] if j is not None else None
        kind = wind_kind(wd, facing)
        tide, nxt = g("sea_level_height_msl"), (m["sea_level_height_msl"][i + 1] if i + 1 < len(m["time"]) else None)
        rows.append({"t": t, "h": tot, "a": round(a, 2), "per": top[1], "wdir": top[2], "sst": g("sea_surface_temperature"),
                     "tide": tide, "tideNext": nxt, "spd": spd, "dir": wd, "kind": kind, "score": score(a, top[1] or 0, Q.get(t), kind, spd, "short", wd, facing), "score_short": score(a, top[1] or 0, Q.get(t), kind, spd, "short", wd, facing), "score_long": score(a, top[1] or 0, Q.get(t), kind, spd, "long", wd, facing), "q": round(Q[t], 2) if Q.get(t) is not None else None,
                     "tmp": x["temperature_2m"][xi[t]] if x and t in xi else None, "wx": x["weather_code"][xi[t]] if x and t in xi else None,
                     "wv": [g("wave_height"), g("wave_period"), g("wave_direction")], "ww": [g("wind_wave_height"), g("wind_wave_period"), g("wind_wave_direction")],
                     "sw": [g("swell_wave_height"), g("swell_wave_period"), g("swell_wave_direction")], "sw2": [g("secondary_swell_wave_height"), g("secondary_swell_wave_period"), g("secondary_swell_wave_direction")],
                     "comps": [[round(h, 2), round(T, 1) if T else None, round(d) if d is not None else None, round(sh[0], 2)] for (h, T, d), sh in zip(comps, shore)]})
    sources = {
        "query_location": {"lat": p["lat"], "lon": p["lon"]},
        "marine": {"provider": "Open-Meteo Marine", "grid": source_grid(marine_response),
                   "parameters": "波浪・うねり・海面水温・海面高度（モデル値）"},
        "tide": {"provider": "Open-Meteo Marine / sea_level_height_msl", "grid": source_grid(marine_response),
                 "description": "平均海面基準のモデル海面高度。潮汐観測所の実測潮位ではありません"},
        "wind": {"provider": wind_model, "grid": source_grid(wind_response),
                 "parameters": "地上10m風速・風向"},
        "weather": {"provider": "Open-Meteo Forecast / 自動選択モデル" if weather_response else "取得なし",
                    "grid": source_grid(weather_response) if weather_response else None,
                    "parameters": "地上2m気温・天気コード"},
    }
    return rows, sources

def main():
    points = json.loads((ROOT / "points.json").read_text(encoding="utf-8"))
    out = {"generated_at": datetime.datetime.now(JST).isoformat(timespec="minutes"), "points": []}
    for p in points:
        try:
            rows, sources = forecast(p)
            out["points"].append({"id": p["id"], "name": p["name"], "lat": p["lat"], "lon": p["lon"],
                                  "facing": p["facing"], "map_url": p.get("map_url"), "facing_source": p.get("facing_source"), "rows": rows, "sources": sources})
        except Exception as e:
            print("skip", p["id"], e, file=sys.stderr)
        time.sleep(0.5)
    if not out["points"]:
        sys.exit("全ポイントで取得に失敗しました")
    d = ROOT / "docs" / "data"; d.mkdir(parents=True, exist_ok=True)
    (d / "forecast.json").write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")

if __name__ == "__main__":
    main()
