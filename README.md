# Twitter Retrosynthesis Bot
[![DOI](https://zenodo.org/badge/330123879.svg)](https://zenodo.org/badge/latestdoi/330123879)

Source code for [@retrosynthchan](https://twitter.com/retrosynthchan)

## Dependencies
- [AiZynthFinder](https://github.com/MolecularAI/aizynthfinder/)
- [Tweepy](https://www.tweepy.org/)

You need to get the [keys and tokens](https://developer.twitter.com/en/docs/authentication/oauth-1-0a) to access Twitter API and modify `settings.json` to run the bot.

## Installation with Docker
1. Install Docker
2. Build the container image and run the image
```
git clone https://github.com/n-yoshikawa/retrosynthesis_bot.git
cd retrosynthesis_bot
docker build -t retrosynth .
docker run -it retrosynth
```
3. Edit `setting.json` to access Twitter API
```
sudo apt install vim  # Install your favorite editor
vim setting.json  # Edit setting.json with your access token
```
4. Run the bot
```
conda activate aizynth-dev
python random_analysis.py  # Run retrosynthesis for a random molecule from PubChem for testing
python bot.py  # Run interactive retrosynthesis bot
```
