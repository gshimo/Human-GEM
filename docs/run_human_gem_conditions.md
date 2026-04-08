# `run_human_gem.py` で外部環境や遺伝子条件を変えながら実行する

このドキュメントは、[`code/examples/run_human_gem.py`](../code/examples/run_human_gem.py) を使って
Human-GEM を実行するときに、外部環境や遺伝子条件をどう変えるかを整理したメモです。

## 1. まず前提

現在の `run_human_gem.py` は、基本的に次の流れだけを行います。

1. `model/Human-GEM.yml` を読み込む
2. 必要なら目的関数を変更する
3. `model.optimize()` を実行する
4. 上位フラックス反応を表示する

そのため、条件を変えたい場合は `model.optimize()` の前に制約を追加します。

```python
model = cobra.io.load_yaml_model(str(model_path))

# ここで環境条件や遺伝子条件を変更する

solution = model.optimize()
```

## 2. 外部環境を変える方法

FBA では、培地や酸素条件は主に exchange 反応の制約として与えます。

### 方法1: `model.medium` を使う

```python
with model:
    medium = model.medium.copy()
    medium["MAR09034"] = 5.0  # glucose uptake を制限
    model.medium = medium

    solution = model.optimize()
    print(solution.status, float(solution.objective_value))
```

`with model:` を使うと、そのブロックを抜けたあとに変更が元に戻るので、
複数条件を試すときに安全です。

### 方法2: exchange 反応の bound を直接変更する

```python
with model:
    rxn = model.reactions.get_by_id("MAR09034")
    rxn.lower_bound = -5.0
    rxn.upper_bound = 1000.0

    solution = model.optimize()
```

一般に exchange 反応では、取り込みは負フラックス、分泌は正フラックスです。
そのため uptake を許可するには `lower_bound` を負に設定します。

### exchange 反応 ID の探し方

Human-GEM では exchange 反応 ID が `EX_glc__D_e` のような名前ではなく、
`MAR...` になっていることがあります。例えば glucose は次のように見つかります。

```bash
uv run python -c "import cobra; m=cobra.io.load_yaml_model('model/Human-GEM.yml'); print([(r.id,r.name) for r in m.exchanges if 'glucose' in (r.name or '').lower()])"
```

このリポジトリのモデルでは、例えば次が確認できます。

- `MAR09034`: `Exchange of glucose`

## 3. 遺伝子条件を変える方法

### 3.1 遺伝子ノックアウト

遺伝子欠損を試すだけなら `knock_out()` が簡単です。

```python
with model:
    model.genes.get_by_id("ENSG00000000419").knock_out()
    solution = model.optimize()
```

遺伝子 ID は Ensembl ID です。

### 3.2 発現低下を近似する

通常の FBA は、発現量そのものを直接読むわけではありません。
発現量を使いたい場合は、それを反応制約に変換する必要があります。

単純な近似として、対象遺伝子に関係する反応の上限と下限を縮める方法があります。

```python
with model:
    gene = model.genes.get_by_id("ENSG00000000419")

    for reaction in gene.reactions:
        if reaction.lower_bound < 0:
            reaction.lower_bound *= 0.3
        if reaction.upper_bound > 0:
            reaction.upper_bound *= 0.3

    solution = model.optimize()
```

これは「発現が低いので通せるフラックスも下げる」という簡易近似です。
ただし、GPR の `and` / `or` を厳密に扱う方法ではありません。

## 4. 複数条件をまとめて回す方法

環境条件と遺伝子条件を同時に変えたい場合は、条件ごとに `with model:` を使って回すと扱いやすくなります。

```python
conditions = [
    {"name": "high_glucose", "glucose": 10.0, "ko": []},
    {"name": "low_glucose", "glucose": 1.0, "ko": ["ENSG00000000419"]},
]

for cond in conditions:
    with model:
        medium = model.medium.copy()
        medium["MAR09034"] = cond["glucose"]
        model.medium = medium

        for gene_id in cond["ko"]:
            model.genes.get_by_id(gene_id).knock_out()

        solution = model.optimize()
        print(cond["name"], solution.status, float(solution.objective_value))
```

## 5. 「遺伝子発現データ」を本格的に使いたいとき

`run_human_gem.py` は単純な FBA 実行例なので、RNA-seq や HPA 由来の発現データをそのまま反映する仕組みは入っていません。

本格的に発現データから組織・細胞型依存のモデルを作りたい場合は、
`code/tINIT/getINITModel2.m` を使う流れが本筋です。
この関数は `arrayData` や `hpaData` を受け取り、発現データに基づいた条件付きモデル構築を行います。

`getINITModel2.m` では、少なくとも次のような入力構造が想定されています。

- `arrayData.genes`
- `arrayData.tissues`
- `arrayData.celltypes`
- `arrayData.levels`
- `arrayData.threshold`

つまり、

- 単純な条件変更: `run_human_gem.py` 側で exchange 制約やノックアウトを入れる
- 発現データベースを使った組織特異的モデル化: `tINIT` を使う

という使い分けになります。

## 6. 実務上のおすすめ

最初は次の順で条件を追加すると扱いやすいです。

1. 目的関数を確認する
2. glucose や amino acids などの exchange 反応を制御する
3. 酸素条件を変える
4. ノックアウトしたい遺伝子を追加する
5. 必要なら複数条件をループで比較する

`run_human_gem.py` 自体に CLI オプションを増やすなら、次のような引数が実用的です。

- `--medium-csv`
- `--set-bound RXN_ID,LB,UB`
- `--knockout-gene GENE_ID`
- `--scale-gene GENE_ID,FACTOR`
