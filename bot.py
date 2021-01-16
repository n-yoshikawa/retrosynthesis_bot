import os
import sys
import glob
import traceback
import subprocess
import urllib.parse
import pandas as pd
from more_itertools import chunked

import tweepy
from rdkit import Chem
from rdkit.Chem import Draw
from aizynthfinder.analysis import ReactionTree

def run(dirname):
    # Run AiZynthFinder CLI
    proc = subprocess.run("aizynthcli --config /home/ubuntu/network/config.yml --output {}/output.hdf5 --smiles {}/smiles.txt > /dev/null".format(dirname, dirname), shell=True, text=True)
    # Generate synthesis route images
    # https://molecularai.github.io/aizynthfinder/html/cli.html#analysing-output
    data = pd.read_hdf("{}/output.hdf5".format(dirname), "table")
    all_trees = data.trees.values
    trees_for_first_target = all_trees[0]
    for itree, tree in enumerate(trees_for_first_target):
        imagefile = f"{dirname}/route{itree:03d}.png"
        ReactionTree.from_dict(tree).to_image().save(imagefile)
    # Combine 4 images into a single image by ImageMagick
    images = glob.glob(f'{dirname}/route*.png')
    for i, imgs in enumerate(chunked(sorted(images), 3)):
        arg = ' '.join(imgs)
        proc = subprocess.run("convert -append -gravity south -splice 0x40 {} {}/result{}.png > /dev/null".format(arg, dirname, i), shell=True, text=True)

class MyStreamListener(tweepy.StreamListener):
    # This function is called every time a tweet containing '@retrosynthchan' is found.
    def on_status(self, status):
        # Avoid reacting retweets
        if hasattr(status, "retweeted_status"):
            return
        # Avoid self reply
        elif status.author.screen_name == '@retrosynthchan':
            return
        else:
            # Get full text (https://docs.tweepy.org/en/latest/extended_tweets.html)
            try:
                text = status.extended_tweet["full_text"]
            except AttributeError:
                text = status.text
        try:
            # SMILES should come right after the screen name
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
        dirname = urllib.parse.quote(smiles, safe='')
        # if dir does not exist (i.e. first run), run retrosynthesis
        if not os.path.exists(dirname):
            try:
                os.mkdir(dirname)
                # Save SMILES as a file to get the result of aizynthfindercli in HDF5 format 
                with open('{}/smiles.txt'.format(dirname), mode='w') as f:
                    f.write(smiles)
                # Generate an image of molecule to attach
                Draw.MolToFile(m, '{}/mol.png'.format(dirname))
                # Trim space by ImageMagick
                proc = subprocess.run("mogrify -trim +repage {}/mol.png > /dev/null".format(dirname), shell=True, text=True)
                # Send a request confirmation message
                res = api.media_upload('{}/mol.png'.format(dirname))
                api.update_status(status='@{} Started retrosynthesis for {}. Please wait for a while.'.format(status.author.screen_name, smiles_org), media_ids=[res.media_id], in_reply_to_status_id=status.id)
                # Run AiZynthFinder
                run(dirname)
            except:
                # If an error occurs, return error message
                print(traceback.format_exc())
                api.update_status(status='@{} Retrosynthesis failed.'.format(status.author.screen_name), in_reply_to_status_id=status.id)

        # Return the result of retrosynthesis
        # Upload up to 4 images
        filenames = glob.glob(f'{dirname}/result*.png')
        media_ids = []
        for filename in sorted(filenames)[:4]:
             res = api.media_upload(filename)
             media_ids.append(res.media_id)
        # If synthesis route images exist, send the result images
        if len(media_ids) != 0:
            api.update_status(status='@{} Retrosynthesis results for {}'.format(status.author.screen_name, smiles_org), media_ids=media_ids, in_reply_to_status_id=status.id)
        # If there is no synthesis route image, send an error message
        else:
            api.update_status(status='@{} No retrosynthesis route was found.'.format(status.author.screen_name), in_reply_to_status_id=status.id)

if __name__ == "__main__":
    # Set Twitter API Access Tokens
    consumer_key = "XXXXXXXXXXXXXXXXXXXXXXXXX"
    consumer_secret = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    access_token = "XXXXXXXXXXXXXXXXXXX-XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    access_secret = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"

    auth = tweepy.OAuthHandler(consumer_key, consumer_secret)
    auth.set_access_token(access_token, access_secret)
    api = tweepy.API(auth)

    # Start stream (https://docs.tweepy.org/en/latest/streaming_how_to.html)
    myStreamListener = MyStreamListener()
    myStream = tweepy.Stream(auth=api.auth, listener=myStreamListener)
    myStream.filter(track=['@retrosynthchan'])
