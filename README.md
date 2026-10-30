# 海点 UMITEN — GitHub Pages 公開用

アップロードされた `umiten-site.zip` をもとに、予報更新とGitHub Pages公開を自動で行うようにしたパッケージです。

## はじめての公開（Windows / ブラウザのみ）

1. [GitHub](https://github.com/) に無料登録・ログイン
2. [New repository](https://github.com/new) から名前 `umiten-site` で `Public` のリポジトリを作成（READMEは追加しない）
3. このZIPをWindowsで「すべて展開」
4. GitHubの空のリポジトリの `uploading an existing file` をクリック（または `Add file` → `Upload files`）
5. 展開したフォルダー **の中身**（`.github`、`docs`、`scripts`、`points.json`、`README.md`）を全選択してドラッグ＆ドロップ
   - フォルダーを丸ごとドラッグすると1階層深くなってしまうので注意
   - 画面下の `Commit changes` で確定
6. リポジトリの `Settings` → 左メニュー `Pages` → `Build and deployment` の `Source` を **GitHub Actions** にする
7. `Actions` → 左側 `update-forecast` → `Run workflow` → `Run workflow` で実行
8. 緑色のチェック（build / deploy）を確認。`Settings` → `Pages` の `Visit site` から公開URLを開く

公開URLの形式: `https://<GitHubユーザー名>.github.io/umiten-site/`

## 中身

- `docs/index.html`：表示画面（HTML/CSS/JavaScript）
- `scripts/forecast.py`：Open-MeteoのAPIからデータ取得・簡易波浪計算
- `points.json`：ポイントの座標と補正パラメータ（現在は3地点で一部仮値）
- `.github/workflows/update.yml`：予報生成と公開を、初回・変更時・3時間ごと・手動で実施
- `docs/data/.gitkeep`：GitHubで空のフォルダーを保持するためのファイル。実データ `forecast.json` は実行時生成

**大切な注意**

- `docs/data/forecast.json` はZIPにありません。Actions成功時に生成され、公開先に含まれます。
- GitHub Pagesの `Source` は **Deploy from a branch** ではなく **GitHub Actions** にしてください。
- GitHub Actionsの設定で `Read and write permissions` に切り替える必要はありません（この版はリポジトリに予報ファイルを書き込みません）。
- 定期実行は最大3時間程度間隔で、GitHub側の都合で遅延することがあります。パブリックリポジトリが60日間操作されないとスケジュールが無効化される場合があります。
- 予報は独自の簡易計算であり、実測により検証されていません。沿岸の潮位や波高が十分正確とは限りません。
- データ提供元 Open-Meteo の無料APIは非商用利用向けで、クレジット表示が必要です。サイトの収益化・商用化前に利用条件を見直してください。
- GitHub Pagesの無料公開には `Public` リポジトリが必要です。中のソースコードも一般公開されます。秘密情報を追加しないでください。
