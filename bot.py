import argparse
import glob
import gc
import os
import json
import subprocess
import traceback
import tweepy
import pandas as pd
import urllib.request
import shutil

from more_itertools import chunked
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

# Set up AiZynthFinder
parser = argparse.ArgumentParser('Retrosynthesis Bot')
parser.add_argument('--config', default='/home/retro/aizynthfinder/model/config.yml')
parser.add_argument('--settings', default='/home/retro/retrosynthesis_bot/settings.json')
args = parser.parse_args()

# Set Twitter API Access Tokens
with open(args.settings) as f:
    d = json.load(f)
    consumer_key = d['consumer_key']
    consumer_secret = d['consumer_secret']
    access_token = d['access_token']
    access_secret = d['access_secret']

auth = tweepy.OAuthHandler(consumer_key, consumer_secret)
auth.set_access_token(access_token, access_secret)
api = tweepy.API(auth)
    
class MyStreamListener(tweepy.StreamListener):
    # This function is called every time a tweet containing '@retrosynthchan' is found.
    def on_status(self, status):
        # Avoid reacting retweets
        if hasattr(status, "retweeted_status"):
            return
        # Avoid self reply
        elif status.author.screen_name == 'retrosynthchan':
            return
        else:
            # Get full text (https://docs.tweepy.org/en/latest/extended_tweets.html)
            try:
                text = status.extended_tweet["full_text"]
            except AttributeError:
                text = status.text
            print(f"@{status.author.screen_name}: {text}")
        if 'media' in status.entities:
            try:
                media = status.extended_entities['media'][0]
                media_url = media['media_url_https']
                filename = media_url.split('/')[-1]
                filename = os.path.join('/tmp', filename)
                urllib.request.urlretrieve(media_url, filename)
                path = os.path.dirname(os.path.abspath(__file__))
                proc = subprocess.run(f"python {path}/get_smiles.py --image {filename}", shell=True, text=True, capture_output=True)
                smiles_org = proc.stdout.rstrip()
            except:
                api.update_status(status='@{} Error: failed to convert the input image into SMILES.'.format(status.author.screen_name), in_reply_to_status_id=status.id)
                return
        else:
            try:
                # SMILES should come right after the screen name
                reply_to = text.split(' ')[0]
                if reply_to != '@retrosynthchan':
                    return
                smiles_org = text.split(' ')[1]
            except:
                # If the input does not contain any spaces, return error message
                api.update_status(status='@{} Error: the input string is invalid.'.format(status.author.screen_name), in_reply_to_status_id=status.id)
                return
        # Check the validity of input SMILES
        m = Chem.MolFromSmiles(smiles_org)
        if m is None:
            # if the input is invalid, return error message
            api.update_status(status='@{} Error: the input string is invalid as SMILES'.format(status.author.screen_name), in_reply_to_status_id=status.id)
            return

        # Get canonical SMILES
        smiles = Chem.MolToSmiles(m)
        # Convert SMILES into safe string (e.g. '/' should not be used as a dir name)
        dirname = Chem.inchi.MolToInchiKey(m)
        # if dir does not exist (i.e. first run), run retrosynthesis
        if os.path.exists(dirname):
            shutil.rmtree(dirname)
        try:
            os.mkdir(dirname)
            # Generate an image of molecule to attach
            drawer = rdMolDraw2D.MolDraw2DCairo(800, 450)
            drawer.DrawMolecule(m)
            drawer.FinishDrawing()
            image_filename = f'{dirname}/mol.png'
            drawer.WriteDrawingText(image_filename)
            # Send a request confirmation message
            res = api.media_upload(image_filename)
            api.update_status(status='@{} Started retrosynthesis for {}. Please wait for a while.'.format(status.author.screen_name, smiles_org), media_ids=[res.media_id], in_reply_to_status_id=status.id)
            # Run AiZynthFinder
            path = os.path.dirname(os.path.abspath(__file__))
            proc = subprocess.run(f"python {path}/generate_images.py --smiles \"{smiles}\" --out_dir {dirname} >> /tmp/aizynthfinder_log.txt 2>&1", shell=True, text=True)
        except:
            # If an error occurs, return error message
            print(traceback.format_exc())
            api.update_status(status='@{} Retrosynthesis failed.'.format(status.author.screen_name), in_reply_to_status_id=status.id)
            return

        # Return the result of retrosynthesis
        # Upload up to 4 images
        filenames = glob.glob(f'{dirname}/result*.png')
        media_ids = []
        for filename in sorted(filenames)[:4]:
             res = api.media_upload(filename)
             media_ids.append(res.media_id)
        # If synthesis route images exist, send the result images
        if len(media_ids) != 0:
            msg = '@{} Retrosynthesis results for '.format(status.author.screen_name)
            tag = ' #ChemicalSpace_Retro'
            if len(smiles_org) > 140-len(msg)-len(tag):
                smiles_trancated = smiles_org[:140-len(msg)-len(tag)-3] + '...'
            else:
                smiles_trancated = smiles_org
            msg += smiles_trancated + tag
            api.update_status(status=msg, media_ids=media_ids, in_reply_to_status_id=status.id)
        # If there is no synthesis route image, send an error message
        else:
            api.update_status(status='@{} No retrosynthesis route was found.'.format(status.author.screen_name), in_reply_to_status_id=status.id)


# Start stream (https://docs.tweepy.org/en/latest/streaming_how_to.html)
print("starting retrosynthesis bot")
myStreamListener = MyStreamListener()
myStream = tweepy.Stream(auth=api.auth, listener=myStreamListener)
myStream.filter(track=['@retrosynthchan'])
