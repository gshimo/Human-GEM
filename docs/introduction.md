# Human-GEM 入門（MATLAB 初心者向け）

このドキュメントは、Human-GEM を MATLAB で使い始めるための最短ガイドです。  
ここまでの会話内容を、実行時に迷いやすい点に絞って整理しています。

## 1. Human-GEM はどんなモデルか

Human-GEM の標準的な使い方は **FBA（Flux Balance Analysis）** です。  
これは「濃度の時間変化」を直接解くモデルではなく、**定常状態のフラックス**を解くモデルです。

- 基本方程式: `S * v = 0`
- `S`: 化学量論行列
- `v`: 反応フラックス

意味:

- 各代謝物について「生成量 - 消費量 = 0（蓄積しない）」を仮定
- その制約の中で、目的関数（例: バイオマス反応）を最大化する

## 2. 「初期値」を与えるモデルではない

通常の FBA で与える入力は、濃度初期値ではなく **制約条件** です。

### 入力（与えるもの）

- 反応の上下限 `lb/ub`
- 交換反応の制約（培地条件、取り込み上限）
- 目的関数 `c`（どの反応を最適化するか）
- 必要に応じてノックアウト制約（遺伝子・反応）

### 出力（得られるもの）

- 反応フラックス `v`
- 目的関数値（例: 成長速度の予測）
- 取り込み・分泌フラックス
- 条件依存の必須遺伝子/必須反応予測

### FBA単体で得られないもの

- 代謝物濃度
- 濃度の時間変化
- 動力学パラメータ（速度定数）ベースの挙動

## 3. MATLAB での最小実行例

```matlab
repo = '/Users/gs/projects/github.com/gshimo/Human-GEM';

% 初回のみ
cd(fullfile(repo,'code'));
HumanGEMInstaller.install;

% モデル読み込み（main ブランチなら .mat 推奨）
matFile = fullfile(repo,'model','Human-GEM.mat');
if exist(matFile,'file')
    tmp = load(matFile);
    fn  = fieldnames(tmp);
    model = tmp.(fn{1});
else
    model = importYaml(fullfile(repo,'model','Human-GEM.yml'), true);
end

% 目的関数設定（Generic human cell biomass reaction）
model = setParam(model,'obj','MAR13082',1);
sol = solveLP(model,1);  % 1 = maximize

fprintf('Objective = %.6f\n', sol.f);
```

## 4. 濃度変化をプロットしたいとき

FBAは定常状態なので、濃度時系列は直接出ません。  
濃度を見たい場合は、交換反応フラックスを外側で時間積分する **簡易 dFBA** を使います。

考え方（外液濃度 `C`）:

- `C(t+dt) = C(t) + v_ex(t) * X(t) * dt`
- `X`: バイオマス濃度
- `v_ex`: 交換反応フラックス（単位: mmol/gDW/h）

このモデルでよく使う交換反応 ID:

- グルコース交換: `MAR09034` (`Exchange of glucose`)
- L-乳酸交換: `MAR09135` (`Exchange of L-lactate`)
- 酸素交換: `MAR09048` (`Exchange of O2`)
- バイオマス反応: `MAR13082`

注意:

- 交換反応は通常「取り込み = 負フラックス」「分泌 = 正フラックス」
- 単位整合（`mmol/L`, `gDW/L`, `h`）を必ず確認する
- より厳密な解析には、専用の dFBA 実装や ODE モデルが必要

## 5. まず何を設定すればよいか（実務上の最小セット）

1. 目的関数（通常は `MAR13082`）
2. 培地条件（主要交換反応の `lb`）
3. 酸素条件（好気/嫌気）
4. 必要なら遺伝子欠損条件

この4点を決めると、最初の FBA 結果（成長・取り込み・分泌）が再現しやすくなります。
