#!/usr/bin/env python3
"""海点 UMITEN の予報データ生成プログラム（Python標準ライブラリのみ）。

処理の流れ
----------
1. Open-Meteo Marine から、沖の波高・周期・波向・うねり・風波などを取得。
2. Open-Meteo Weather から、10 m 高度の風速・風向と直近24時間の予報履歴を取得。
3. 風向と各海岸の向きを比較し、その時刻の「面が整いやすい条件」を算出。
4. 過去の風により面が荒れていた記憶を指数緩和モデルで引き継ぎ、
   1時間ごとに面の状態 Q（0〜1）を更新。
5. Q を板種別総合スコアに入力し、docs/data/forecast.json を出力。
6. point別の「scale」係数で、近隣ポイントの砕波高の系統差を暫定補正。
   海の子はつるかめより少しサイズが大きいという現地情報を参考に調整。

【重要】面モデルは物理法則そのものではなく、予報用の経験的モデルです。
海岸近くの風・砕波や潮による地形変化を直接観測していないため、
現地観測や波の映像と突き合わせて係数を校正する必要があります。
"""
import json, math, os, ssl, sys, time, datetime, pathlib, urllib.parse, urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
G = 9.81
JST = datetime.timezone(datetime.timedelta(hours=9))
# 面の指数緩和時間の基準値 [hour]。後で風向・風速により可変にする。
# Q が目標値に63.2%近づく時間が「時定数」であり、完全に整う時間ではない。
LAG_UP = 3.0       # 面が整う（Q上昇）ときの基準時定数
LAG_DOWN = 1.8     # 面が荒れる（Q低下）ときの基準時定数
# 強いオンショアが続けば荒れやすく、適度なオフショアなら早く回復する。
# 一方、急な風向変化でQを即座に目標値へ飛ばさないことが狙い。
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

def clamp(value, low, high):
    """値を指定区間に収める。特に面Qは必ず 0〜1 に抑える。"""
    return max(low, min(high, value))


def wind_components(speed, direction_from, facing):
    """10 m 風速を岸に直交する成分・海岸に平行な成分に分解する。

    Open-Meteo の wind_direction_10m は「風が吹いてくる方角」。
    facing は「岸から海を見たときの向き」（例: 海が東なら90度）。

    d = 0度: 海から岸へのオンショア、d = 180度: 陸から海へのオフショア。
    onshore  = max(V*cos(d), 0)
    offshore = max(-V*cos(d), 0)
    cross    = V*abs(sin(d))

    方角を5段階に区切るのではなく、連続値を使うことで、
    例えば風向が59度→61度と変わっても突然評価が飛ばないようにする。
    """
    if speed is None or facing is None:
        return None
    speed = max(0.0, float(speed))
    # 0m/sでは風向は意味を持たない。APIが風向をnullで返しても、
    # 無風として面の自然回復を計算できるようにする。
    if direction_from is None:
        return {"on": 0.0, "off": 0.0, "side": 0.0, "speed": speed} if speed <= 0.5 else None
    d = math.radians(float(direction_from) - float(facing))
    c, s = math.cos(d), math.sin(d)
    return {"on": speed * max(c, 0.0),
            "off": speed * max(-c, 0.0),
            "side": speed * abs(s), "speed": speed}


def wind_sea_fraction(marine_hour):
    """全波に含まれる風波のエネルギー比を [0,1] で近似する。

    波のエネルギーは波高の二乗に概ね比例するため、H_wind^2/H_total^2を使う。
    海面そのものの波状の細かさを直接測った値ではなく、補助情報にすぎない。
    モデル分解能が粗く、風波の値が不安定な場合もあるので、寄与は小さくする。
    """
    if not marine_hour:
        return None
    total, wind = marine_hour.get("wave_height"), marine_hour.get("wind_wave_height")
    if total is None or wind is None or total <= 0.05 or wind < 0:
        return None
    return clamp((wind / total) ** 2, 0.0, 1.0)


def surface_equilibrium(components, gust=None, marine_hour=None):
    """仮に同じ風が長時間続いたときの「面の落ち着き先」Q_eqを見積もる。

    Q_eqが高いほど面が整いやすい。計算は滑らかな飽和関数を使う。

      Q_eq = 0.96
             - 0.86 [1-exp(-(V_on/3.8)^1.55)]
             - 0.50 [1-exp(-(V_side/4.5)^1.6)]
             - 強風・ガスト・風波の小補正

    風速や方向による寄与は実測に基づく確定係数ではなく、海点の暫定値。
    特に「うねり周期が長いから面は必ず綺麗」という誤った仮定はしない。

    ここで決まるのは「風が十分長く続いた先の状態」であり、
    現時刻の面Qは別の surface_quality で履歴を考慮して求める。
    """
    on, side, speed = components["on"], components["side"], components["speed"]
    on_penalty = 0.86 * (1.0 - math.exp(-((on / 3.8) ** 1.55)))
    side_penalty = 0.50 * (1.0 - math.exp(-((side / 4.5) ** 1.6)))
    # 8m/s以上の風は、純粋なオフショアでも飛沫や風による影響が目立つ。
    # ただし「波が抑えられてテイクオフしにくい」のは総合評価側でも考慮する。
    strong_penalty = 0.12 * (1.0 - math.exp(-((max(0.0, speed - 7.0) / 3.0) ** 1.5)))

    # 強い瞬間風速の影響は、オン・サイド成分があるときに面へ反映。
    # 最大0.05点(0〜1スケール)なので、ガスト予報の不確実性に過度に依存しない。
    gust_penalty = 0.0
    if gust is not None and gust >= speed:
        excess = max(0.0, float(gust) - speed - 2.0)
        gust_penalty = 0.05 * clamp(excess / 5.0, 0.0, 1.0) * clamp(on / 3.0 + side / 6.0, 0.0, 1.0)

    fraction = wind_sea_fraction(marine_hour)
    # 比率が高いときは乱れが残りやすい。ただし寄与は最大0.08に抑える。
    sea_penalty = 0.08 * fraction if fraction is not None else 0.0
    return clamp(0.96 - on_penalty - side_penalty - strong_penalty
                 - gust_penalty - sea_penalty, 0.04, 0.96)


def surface_time_constant(previous_q, target_q, components, marine_hour,
                          up=LAG_UP, down=LAG_DOWN):
    """現在の面Qが目標Q_eqに近づく時定数 tau[hour] を決める。

    ・荒れる場合: オンショアが強いほど速く荒れる（tau小）。
    ・整う場合: 無風でも徐々に回復。適度なオフショアなら少し速い。
    ・長周期うねり: 回復を最大15%だけ促進する弱い補助項にする。
    ・風波が多い: 残ったチョッピーな波を考え、回復を最大25%遅らせる。

    元案の「荒れの蓄積度/(オフショア風速×周期)」は使わない。
    ゼロ風速で発散する問題と、長周期の影響を過大評価する問題を避ける。
    """
    on, off, side = components["on"], components["off"], components["side"]
    if target_q < previous_q:
        # 風が横から当たるだけでも荒れるが、オンショアほどは速くない。
        tau = down / (1.0 + 0.22 * on + 0.10 * side)
        return clamp(tau, 0.60, max(0.60, float(down)))

    # オフショアは最大6m/sで飽和させる（強いほど速く綺麗、とはしない）。
    tau = up / (1.0 + 0.16 * min(6.0, off))
    swell_period = marine_hour.get("swell_wave_period") if marine_hour else None
    if swell_period is not None and swell_period > 0:
        # 8秒以下なら補正なし、14秒以上でも効果は15%まで。
        tau *= 1.0 - 0.15 * clamp((swell_period - 8.0) / 6.0, 0.0, 1.0)
    fraction = wind_sea_fraction(marine_hour)
    if fraction is not None:
        tau *= 1.0 + 0.25 * fraction
    return clamp(tau, 1.0, max(1.0, 1.5 * float(up)))


def surface_quality(wind_times, speeds, directions, facing,
                    gusts=None, marine_by_time=None, up=LAG_UP, down=LAG_DOWN):
    """風の履歴から「面の整い具合 Q」を逐次計算する（0〜1）。

    原理：各時刻の目標 Q_eq と、前時刻の面 Q_old を使い、

      Q_new = Q_eq + (Q_old - Q_eq) * exp(-Δt/tau)

    で更新する（指数緩和、一次遅れ）。
    Δt=tau なら、目標値との差が約37%まで縮む（63%の変化）。
    直前までオンショアで面が荒れていた場合、オフショアになっても
    すぐQ_eqまで跳ね上がらず、数時間かけて自然に回復する。

    入力には過去1日の予報アーカイブを含める。予報履歴≠現地観測であり、
    風向が急変する沿岸では現実の状態とずれる可能性がある。
    欠測した時刻は面Qを据え置き、無理に「完全に整った」とは判断しない。

    返り値: [(Q, Q_eq, tau), ...]。Qだけが公開画面の面評価に用いられる。
    Q_eq/tauは計算の妥当性を後から調査するためJSONにも保存する。
    """
    gusts = gusts if gusts is not None else [None] * len(wind_times)
    marine_by_time = marine_by_time or {}
    result = []
    q_previous = None
    previous_timestamp = None

    for i, iso_time in enumerate(wind_times):
        parts = wind_components(speeds[i], directions[i], facing)
        stamp = datetime.datetime.fromisoformat(iso_time)
        dt = (stamp - previous_timestamp).total_seconds() / 3600 if previous_timestamp else 1.0
        previous_timestamp = stamp
        marine_hour = marine_by_time.get(iso_time)

        if parts is None:
            # 未取得は「オンショア」とみなさない。前のQを保持し、診断値は欠測。
            result.append((q_previous, None, None))
            continue

        gust = gusts[i] if i < len(gusts) else None
        target = surface_equilibrium(parts, gust, marine_hour)
        if q_previous is None:
            # 計算系列の先頭。過去の風が不明なので最初は目標値を仮初期値にする。
            # 過去24時間分を先読みしているため、今日の表示時点では初期値の影響が小さい。
            q_previous = target
            tau = None
        else:
            tau = surface_time_constant(q_previous, target, parts, marine_hour, up, down)
            # Open-Meteo hourlyは通常1時間刻み。異常な日時や重複を防御。
            dt = clamp(dt, 0.0, 6.0)
            q_previous = target + (q_previous - target) * math.exp(-dt / tau)
            q_previous = clamp(q_previous, 0.0, 1.0)
        result.append((q_previous, target, tau))
    return result


def wind_score_cap(kind, spd):
    """現在の風向・風速から、その時間帯のスコア上限を設定する。

    面の良さ Q は過去の風を考慮するが、現在オンショアに変わった場合に
    古い良好な Q だけで高得点にならないようにする。試験的な経験則。
    """
    if spd is None or not math.isfinite(float(spd)):
        return None
    # 風向が欠測していても、0.5m/s以下なら影響が十分小さいので減点不要。
    # 風向欠測の有風時は推測せず欠測扱いとする。
    if spd <= 0.5:
        return 100
    if kind is None:
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


# 波高は「良い／悪いを単独決定する閾値」ではなく、板種に対する適性の一要素。
# 2026-11-02: 小さい波は一律×というハードキャップを撤廃する。
# 波高の補正を段階的に効かせ、実際の点数は面Q（過去の風履歴込み）、
# 風向・風速、周期と合算して決める。特にロングはヒザ〜モモでも
# 面が整い風が穏やかなら、△を十分に狙えるよう適合度を広くする。
# ※波高値とサイズ名称の対応は目安。砕波形状・地形によって異なる。
# ※係数は経験則の仮設定で、ポイントの実測による統計校正は未実施。
SHORT_WAIST_FIT_BOOST = 0.06
SHORT_SIZE_FIT_POINTS = (
    (0.00, 0.00), (0.30, 0.10), (0.50, 0.28),
    (0.70, 0.73), (0.90, 0.96), (1.10, 1.00),
    (1.50, 1.00), (1.70, 0.97), (1.90, 0.87),
    (2.15, 0.56), (2.50, 0.18), (3.00, 0.03),
)
# ロングの下限を広げる：ヒザ付近0.4 mでも適性0.49、
# モモ付近0.5 mなら0.66。波が小さくても「面が整い、周期がある」
# ケースでは、もとのガウス型より不当に低評価にならない。
LONG_SIZE_FIT_POINTS = (
    (0.00, 0.00), (0.20, 0.12), (0.30, 0.30),
    (0.40, 0.49), (0.50, 0.66), (0.65, 0.82),
    (0.85, 0.96), (1.08, 1.00), (1.30, 0.88),
    (1.50, 0.62), (1.70, 0.30), (2.00, 0.10),
    (2.40, 0.02), (3.00, 0.00),
)
BOARD_PROFILES = {
    'short': {'fit_points': SHORT_SIZE_FIT_POINTS, 'size_weight': 0.45},
    'long':  {'fit_points': LONG_SIZE_FIT_POINTS,  'size_weight': 0.30},
}


def smoothstep(x):
    """0〜1のなめらかな補間。区切りの境界で採点が飛ばないようにする。"""
    u = clamp(x, 0.0, 1.0)
    return u * u * (3.0 - 2.0 * u)


def wave_presence_factor(height):
    """ほとんど波がないときだけ減点する、連続的な「波の存在」補正。

    波高がちょうど0 mでも面や風だけで△になることを避けたい一方、
    スネ→ヒザなどの区分ごとに点数を切ってしまうのも不自然。
    0.12〜0.38 mの範囲だけをsmoothstepで0〜1に移行させる。
    0.38 m以上なら1.0であり、ヒザ〜モモのロングを妨げない。
    """
    return smoothstep((max(0.0, height) - 0.12) / 0.26)


def size_suitability(a, board='short'):
    """サイズ適合度 F: ショート／ロングの別に滑らかな曲線で評価する。

    このFは「何点以下」という上限ではない。最終点は風・風履歴・
    面・周期の同時条件で決まる。ロングはヒザ〜モモも遊べる適性として扱う。
    """
    if a is None:
        return None
    height = max(0.0, float(a))
    pts = BOARD_PROFILES[board]['fit_points']
    if height <= pts[0][0]:
        return pts[0][1]
    for (h0, y0), (h1, y1) in zip(pts, pts[1:]):
        if height <= h1:
            weight = smoothstep((height - h0) / (h1 - h0))
            fit = y0 + (y1 - y0) * weight
            if board == 'short' and 0.70 < height < 0.90:
                # 従来のコシ補正は維持する。両端で補正量ゼロになる。
                x = (height - 0.70) / 0.20
                fit += SHORT_WAIST_FIT_BOOST * math.sin(math.pi * x) ** 2
            return clamp(fit, 0.0, 1.0)
    return pts[-1][1]


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


def period_score_points(period, board="short"):
    """周期の寄与点（0〜10点）を返す。総合点の強弱とは別に上限判定がある。

    【変更の背景：2026-10-31 とっぷさんての現地情報】
    ムネ前後・西風の条件でも、約7秒の波が一律3点で評価されると、
    面がそこそこ整っていても▼になりやすかった。
    7秒前後でも遊べる波があるため、ショートに限って6〜10秒の
    周期点を**連続的に**変化させる。

    ショート：
      T < 6 s   : 0点
      6〜8 s    : 3点 → 7点（直線補間。7秒なら5点）
      8〜10 s   : 7点 → 10点（直線補間）
      T >= 10 s : 10点

    ロング：既存の段階点（6秒未満0、6〜8秒3、8〜10秒7、10秒以上10）を維持。

    【重要】周期8秒未満は後段の period_cap=69 を維持するため、
    この補正で○が頻発することはない。7秒でも、面・風・サイズが
    悪ければ▼や×になり得る。観測で確定した係数ではない。
    """
    if board == "long":
        return 10.0 if period >= 10 else 7.0 if period >= 8 else 3.0 if period >= 6 else 0.0
    if period < 6:
        return 0.0
    if period < 8:
        return 3.0 + 2.0 * (period - 6.0)
    if period < 10:
        return 7.0 + 1.5 * (period - 8.0)
    return 10.0


def score(a, per, q, kind=None, spd=None, board="short", wind_dir=None, facing=None):
    """総合点0〜100。サイズだけで等級を決めず、風履歴・面・周期を中心に評価。

    S0 = (70 * Q**2.2 + 20 * F + P) * (1 - w*(1-F)) * 0.92
      Q : 面の整い具合。オンショア履歴、回復時定数、ガスト、風波を反映。
      F : 板種別の波高適合度。サイズの単独閾値やハードキャップは使用しない。
      P : 周期の加点0〜10。
      w : サイズ重み。ショート0.45、ロング0.30。
    S = S0 * offshore_suitability * wave_presence_factor を計算し、
    風向・風速および周期の上限で制限する。

    サイズ重みは旧式0.65より小さくし、サイズ由来の二重三重減点を解消。
    ロングはヒザ〜モモでも、十分なQと好風向・周期なら△になり得る。
    一方、実質波がないときは wave_presence_factor で連続的に減点する。
    ◎・○の境界85/70などは従来通り。頻度の統計校正は未実施。
    """
    # 欠測と「実際に0だった値」を分離する。
    if any(v is None for v in (a, per, q, spd)):
        return None
    if not all(math.isfinite(float(v)) for v in (a, per, q, spd)) or per <= 0:
        return None
    if kind is None and spd > 0.5:
        return None
    a = max(0.0, float(a))
    q = clamp(float(q), 0.0, 1.0)
    fit = size_suitability(a, board)
    period_points = period_score_points(per, board)

    # 新式：サイズは滑らかな適性に過ぎない。サイズのハードキャップはかけない。
    size_weight = BOARD_PROFILES[board]['size_weight']
    raw = (70 * q ** 2.2 + 20 * fit + period_points)
    raw *= (1.0 - size_weight * (1.0 - fit)) * 0.92
    raw *= wave_presence_factor(a)

    # 風の向きは直前までの「面の荒れ」とは独立に、現在の吹き方も評価する。
    alignment = offshore_suitability(wind_dir, facing, spd, board)
    if alignment is not None:
        raw *= alignment
    # 周期の上限は維持。6〜8秒でも遊べるが○を濫発させない。
    period_cap = 49 if per < 6 else 69 if per < 8 else 100
    wind_cap = wind_score_cap(kind, spd)
    if wind_cap is None:
        return None
    return max(0, min(100, round(min(raw, period_cap, wind_cap))))


def source_grid(payload):
    """Open-Meteo API の応答に含まれる代表格子点座標。観測所の位置ではない。"""
    lat, lon = payload.get("latitude"), payload.get("longitude")
    return {"lat": round(float(lat), 5) if lat is not None else None,
            "lon": round(float(lon), 5) if lon is not None else None}


def forecast(p):
    q = urllib.parse.urlencode({"latitude": p["lat"], "longitude": p["lon"], "timezone": "Asia/Tokyo", "forecast_days": 7})
    # 風が変わる前の海況（風波比率）も追えるよう海洋予報も過去1日を取得する。
    # past_days は「過去のモデル予報値」であり、ブイ観測の履歴ではない。
    marine_response = get(f"https://marine-api.open-meteo.com/v1/marine?{q}&hourly={MARINE}&past_days=1")
    m = marine_response["hourly"]
    # ガスト（風の瞬間最大値）も取得できれば面に小さく反映する。
    # JMA Seamlessの変数提供状況は時期・地域で異なり得るため段階的にフォールバック。
    # ①JMA+ガスト → ②JMAの風のみ → ③自動選択+ガスト → ④自動選択の風のみ。
    # ガストが無くても、面のモデル自体は計算可能にする。
    wind_base = f"https://api.open-meteo.com/v1/forecast?{q}&wind_speed_unit=ms&past_days=1"
    choices = [
        ("wind_speed_10m,wind_direction_10m,wind_gusts_10m", "jma_seamless", "JMA Seamless"),
        ("wind_speed_10m,wind_direction_10m", "jma_seamless", "JMA Seamless（ガストなし）"),
        ("wind_speed_10m,wind_direction_10m,wind_gusts_10m", None, "自動選択モデル"),
        ("wind_speed_10m,wind_direction_10m", None, "自動選択モデル（ガストなし）"),
    ]
    last_error = None
    for variables, model_arg, model_label in choices:
        url = f"{wind_base}&hourly={variables}"
        if model_arg:
            url += f"&models={model_arg}"
        try:
            wind_response = get(url)
            w = wind_response["hourly"]
            # 一部のモデルは先の時間帯で欠測することがある。
            # 「最初の1点だけ風速があれば成功」とすると週間の総合点が欠けるため、
            # Marine側の表示予定時刻の90%以上に風向・風速があるか確認する。
            # カバー率が足りなければ、次の候補モデルを試す。
            today_jst = datetime.datetime.now(JST).strftime("%Y-%m-%d")
            future_times = {t for t in m["time"] if t[:10] >= today_jst}
            available_times = {
                t for t, speed, direction in zip(
                    w["time"], w["wind_speed_10m"], w["wind_direction_10m"]
                ) if speed is not None and (direction is not None or speed <= 0.5)
            }
            coverage = len(future_times & available_times) / max(1, len(future_times))
            if coverage < 0.90:
                raise ValueError(f"週間予報で風データが不足: カバー率{coverage:.0%}")
            wind_model = "Open-Meteo / " + model_label
            break
        except Exception as exc:
            last_error = exc
    else:
        raise RuntimeError(f"風の予報取得に失敗: {last_error}")
    wi = {t: i for i, t in enumerate(w["time"])}
    # 天気・気温も今日の0時から表示する。past_days=1で時刻の欠落を防ぐ。
    # これはモデルの過去時刻の予報値であり、観測記録ではない。
    try:   # 天気と気温(表示用。取れなくても続行)
        weather_response = get(f"https://api.open-meteo.com/v1/forecast?{q}&hourly=temperature_2m,weather_code&past_days=1")
        x = weather_response["hourly"]
        xi = {t: i for i, t in enumerate(x["time"])}
    except Exception:
        weather_response, x, xi = None, None, {}
    # 岸の波高補正係数 scale はポイントごとの調整用。
    # 例：大貫海の子はつるかめよりややサイズが出やすいというユーザーの
    # 現地情報を反映し、基準0.75に対して0.79（約5.3%増）とする。
    # これは単独の実測で確立した係数ではなく、地形・波向の違いが
    # 入っている屈折計算と重ね掛けされる。必ず海の子が大きいとは限らない。
    # 今後は同時刻の両地点の波高記録から scale を校正する。
    facing, depth, scale = p["facing"], p.get("depth_m", 4.0), p.get("scale", 1.0)
    up, down = p.get("lag_up_h", LAG_UP), p.get("lag_down_h", LAG_DOWN)
    # MarineとWeatherでは系列の開始時刻や格子間隔が異なることがある。
    # 配列番号で対応させず、必ず現地時刻（ISO時刻文字列）で結合する。
    marine_by_time = {
        t: {key: m[key][i] for key in ("wave_height", "wind_wave_height", "swell_wave_period")}
        for i, t in enumerate(m["time"])
    }
    state = surface_quality(
        w["time"], w["wind_speed_10m"], w["wind_direction_10m"], facing,
        gusts=w.get("wind_gusts_10m"), marine_by_time=marine_by_time,
        up=up, down=down,
    )
    Q = dict(zip(w["time"], state))   # 値は(Q, Q_eq, tau)のタプル
    # 可視化で「面」を表示するときはQ[時刻][0]だけを使う。
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
        # 波高の欠測(None)と実際の0mを混同しない。欠測を0m扱いすると
        # 「データがないだけなのに×評価」が起こるため、ここで区別する。
        tot = g("wave_height")
        rest = math.sqrt(max(0.0, tot * tot - sum(c[0] ** 2 for c in comps))) if tot is not None else 0.0
        if rest > 0.05:
            comps.append((rest, g("wave_period"), g("wave_direction")))
        shore = [(to_shore(h, T, d, facing, depth), T, d) for h, T, d in comps]
        a = (scale * math.sqrt(sum(s[0] ** 2 for s in shore)) if shore else 0.0) if tot is not None else None
        top = max(shore, key=lambda s: s[0]) if shore else (0, g("wave_period"), g("wave_direction"))
        # 主成分の周期が欠測しても全波の周期がある場合はそちらを使う。
        # 両方なければ「0秒の短周期」と偽らず、総合点を欠測にする。
        period = top[1] if top[1] is not None and top[1] > 0 else g("wave_period")
        if period is not None and period <= 0:
            period = None
        j = wi.get(t)
        spd = w["wind_speed_10m"][j] if j is not None else None
        wd = w["wind_direction_10m"][j] if j is not None else None
        kind = wind_kind(wd, facing)
        # 同時刻のQ（履歴反映済み）とQ_eq（今の風が続いた場合の目標値）を分離。
        state_value = Q.get(t, (None, None, None))
        q_actual, q_target, q_tau = state_value
        tide, nxt = g("sea_level_height_msl"), (m["sea_level_height_msl"][i + 1] if i + 1 < len(m["time"]) else None)
        # 採点はショート/ロングそれぞれ1回ずつ。表示箇所に応じて採点を
        # 繰り返すと、後で式が分岐した際に値が不一致になるため避ける。
        short_score = score(a, period, q_actual, kind, spd, "short", wd, facing)
        long_score = score(a, period, q_actual, kind, spd, "long", wd, facing)
        rows.append({"t": t, "h": tot, "a": round(a, 2) if a is not None else None, "per": period, "wdir": top[2], "sst": g("sea_surface_temperature"),
                     "tide": tide, "tideNext": nxt, "spd": spd, "dir": wd, "kind": kind, "score": short_score, "score_short": short_score, "score_long": long_score,
                     "q": round(q_actual, 3) if q_actual is not None else None,
                     # 検証用: HTMLは無視してよいが、時系列でタイムラグを確認できる。
                     "q_target": round(q_target, 3) if q_target is not None else None,
                     "q_tau_h": round(q_tau, 2) if q_tau is not None else None,
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
                 "parameters": "地上10m風速・風向・ガスト（提供される場合）／直近1日のモデル履歴を面評価へ使用"},
        "weather": {"provider": "Open-Meteo Forecast / 自動選択モデル" if weather_response else "取得なし",
                    "grid": source_grid(weather_response) if weather_response else None,
                    "parameters": "地上2m気温・天気コード"},
    }
    sources["surface_model"] = {
        "name": "面の一次遅れモデル v2（経験則、未検証）",
        "description": "風向・風速・ガスト・風波比率から目標面Qを算出。風の履歴と指数緩和時間を使い、面の荒れ・回復の時間差を表現する。",
        "parameters": {"baseline_recovery_h": up, "baseline_roughening_h": down},
        "notice": "風・波のモデル格子値に基づく推定であり、現地の面状態の実測ではない"
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
