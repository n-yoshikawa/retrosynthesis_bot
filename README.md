# Twitter Retrosynthesis Bot
Source code for [@retrosynthchan](https://twitter.com/retrosynthchan)

## Dependencies
- [AiZynthFinder](https://github.com/MolecularAI/aizynthfinder/)
- [Tweepy](https://www.tweepy.org/)

## How to use
1. Install dependencies
2. Get [Access Tokens](https://developer.twitter.com/ja/docs/basics/authentication/guides/access-tokens) and modify `settings.json`
3. Run the bot

```
python --config {path to config.yml of AiZynthFinder} --settings {path to settings.json} bot.py
```

## Set Up Example
The bot is running on AWS t2.large instances with Ubuntu Server 20.04 LTS.

```
# Install AiZynthFinder
sudo apt-get install libxrender1
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
exec "$SHELL"
conda install python=3.7
conda install -c rdkit "rdkit=>2019.09.1" -y
conda install -c anaconda "tensorflow>=2.1.0" -y
conda install -c anaconda graphviz -y
git clone https://github.com/MolecularAI/aizynthfinder.git
cd aizynthfinder
python -m pip install -e .
mkdir model
python aizynthfinder/tools/download_public_data.py model
aizynthcli --config model/config.yml --smiles 'c1ccccc1'
# Install the bot
cd ~
pip install tweepy
pip install pubchempy
git clone https://github.com/n-yoshikawa/retrosynthesis_bot.git
cd retrosynthesis_bot/
vim settings.json  # Set Twitter API Key
python random_analysis.py  # Run retrosynthesis for random PubChem molecule
python bot.py  # Run retrosynthesis bot
```
