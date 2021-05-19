import argparse
import glob
import json
import os
import random

import tweepy
import pubchempy as pcp

from generate_images import generate_images

# Set up AiZynthFinder
parser = argparse.ArgumentParser('Random Retrosynthesis Bot')
parser.add_argument('--config', default='/root/aizynthfinder/model/config.yml')
parser.add_argument('--settings', default='/root/retrosynthesis_bot/settings.json')
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

# Loop until retrosynthesizable molecule is found
while True:
    cid = random.randint(10000, 1000000)
    compound = pcp.get_compounds(cid, 'cid')[0]
    dirname = compound.inchikey
    if not os.path.exists(dirname):
        generate_images(compound.isomeric_smiles, dirname, args.config)
    # Upload up to 4 images
    filenames = glob.glob(f'{dirname}/result*.png')
    media_ids = []
    for filename in sorted(filenames)[:4]:
         res = api.media_upload(filename)
         media_ids.append(res.media_id)
    # If synthesis route images exist, send the result images
    if media_ids:
        if compound.synonyms:
            name = compound.synonyms[0]
        else:
            name = compound.iupac_name
        api.update_status(status=f'Retrosynthetic analysis of PubChem CID {cid}. https://pubchem.ncbi.nlm.nih.gov/compound/{cid}', media_ids=media_ids)
        break
