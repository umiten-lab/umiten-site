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
    return max(q, .6) if spd <= 1.5 else q            # ほぼ無風ならベタ凪に近い

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

def score(a, per, q):
    """風(面)を最重視。サイズ30点、面55点、周期15点。波が小さい日は面と周期の点を割り引く。"""
    if a is None or q is None:
        return None
    pp = 15 if per >= 8 else 8 if per >= 6 else 0
    return round(min(a, 1.5) / 1.5 * 30 + (q * 55 + pp) * min(1, a / 0.4))

def forecast(p):
    q = urllib.parse.urlencode({"latitude": p["lat"], "longitude": p["lon"], "timezone": "Asia/Tokyo", "forecast_days": 7})
    m = get(f"https://marine-api.open-meteo.com/v1/marine?{q}&hourly={MARINE}")["hourly"]
    wq = f"https://api.open-meteo.com/v1/forecast?{q}&hourly=wind_speed_10m,wind_direction_10m&wind_speed_unit=ms&past_days=1"   # 前日分は、風向きの変化の履歴(面のラグ)の計算用
    try:
        w = get(wq + "&models=jma_seamless")["hourly"]          # 気象庁モデル
        assert any(v is not None for v in w["wind_speed_10m"])
    except Exception:
        w = get(wq)["hourly"]
    wi = {t: i for i, t in enumerate(w["time"])}
    try:   # 天気と気温(表示用。取れなくても続行)
        x = get(f"https://api.open-meteo.com/v1/forecast?{q}&hourly=temperature_2m,weather_code")["hourly"]
        xi = {t: i for i, t in enumerate(x["time"])}
    except Exception:
        x, xi = None, {}
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
                     "tide": tide, "tideNext": nxt, "spd": spd, "dir": wd, "kind": kind, "score": score(a, top[1] or 0, Q.get(t)), "q": round(Q[t], 2) if Q.get(t) is not None else None,
                     "tmp": x["temperature_2m"][xi[t]] if x and t in xi else None, "wx": x["weather_code"][xi[t]] if x and t in xi else None,
                     "wv": [g("wave_height"), g("wave_period"), g("wave_direction")], "ww": [g("wind_wave_height"), g("wind_wave_period"), g("wind_wave_direction")],
                     "sw": [g("swell_wave_height"), g("swell_wave_period"), g("swell_wave_direction")], "sw2": [g("secondary_swell_wave_height"), g("secondary_swell_wave_period"), g("secondary_swell_wave_direction")],
                     "comps": [[round(h, 2), round(T, 1) if T else None, round(d) if d is not None else None, round(sh[0], 2)] for (h, T, d), sh in zip(comps, shore)]})
    return rows

def main():
    points = json.loads((ROOT / "points.json").read_text(encoding="utf-8"))
    out = {"generated_at": datetime.datetime.now(JST).isoformat(timespec="minutes"), "points": []}
    for p in points:
        try:
            out["points"].append({"id": p["id"], "name": p["name"], "lat": p["lat"], "lon": p["lon"], "facing": p["facing"], "rows": forecast(p)})
        except Exception as e:
            print("skip", p["id"], e, file=sys.stderr)
        time.sleep(0.5)
    if not out["points"]:
        sys.exit("全ポイントで取得に失敗しました")
    d = ROOT / "docs" / "data"; d.mkdir(parents=True, exist_ok=True)
    (d / "forecast.json").write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")

if __name__ == "__main__":
    main()
