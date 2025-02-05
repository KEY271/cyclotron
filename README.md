# 卓上サイクロトロンシミュレータ

卓上サイズのサイクロトロンのシミュレータです。

## 動かす前に
OpenFOAM をインストールしておいてください。
- OpenFOAM には Foundation 版と ESI 版がありますが、必ず後者をインストールしてください。
- 最新の OpenFOAM にはなぜか cfMesh(メッシュ作成に使う)が付属しないので、少し古いですが v2312 をインストールしてください。

`openscad`ディレクトリの中の`.stl`ファイルは`.scad`ファイルから作られています。これを自分で作る場合は OpenSCAD という3Dモデルを作成するソフトをインストールしてください。

## 動かしてみる
まず`particleSim/physicalProperties.sample`をコピーして`particleSim/physicalProperties`を作ってください。
- これは陽イオンの運動のシミュレーションの設定ファイルです。

次に`openfoam2312`を実行してOpenFOAMのシェルセッションを開始し、以下のシェルスクリプトを順番に実行してください。
- `./Allrun.pre`: メッシュを作成します。
- `./Allrun.main`: `electrostaticFoam`で静電場のシミュレーションを行います。
- `./Allrun.sim`: 陽イオンの運動のシミュレーションを行います。
    - `./Allrun.sim: line 8: 98342 Abort trap: 6`というエラーが出た場合は`source ./Allrun.sim`としてみてください(macだと必要？)。

無事に終了すると`results`ディレクトリに結果が出力されます。結果ファイルは次のような構成になっています。

|項目|説明|
|---|---|
|b|磁場、単位はテスラ|
|count|検出器に当たった粒子の数|
|fallen|フィラメント・導入線に当たった粒子の数|
|r_hit|真空容器・Dee電極に当たった粒子の数(動径方向)|
|z_hit|真空容器・Dee電極に当たった粒子の数(鉛直方向)|

`particleSim/physicalProperties`の設定で`trajectory`を`true`にしていた場合は`results/trajectory`に軌跡が出力されます。

## 他のファイル・ディレクトリの説明
- `./Allclean`: 途中生成ファイルなどを全て削除します。
- `electrostaticFoam_P`: 陽イオンのシミュレーションをするソルバーです。
