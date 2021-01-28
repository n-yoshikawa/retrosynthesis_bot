import argparse
import glob
import json
import os
import random

import tweepy
import pubchempy as pcp

from more_itertools import chunked
from aizynthfinder.aizynthfinder import AiZynthFinder

import utils

# Set up AiZynthFinder
parser = argparse.ArgumentParser('Random Retrosynthesis Bot')
parser.add_argument('--config', default='/home/ubuntu/network/config.yml')
parser.add_argument('--settings', default='/home/ubuntu/settings.json')
args = parser.parse_args()

finder = AiZynthFinder(configfile=args.config)
finder.stock.select("zinc")
finder.expansion_policy.select("uspto")
finder.filter_policy.select("uspto")

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
    print(cid, compound.isomeric_smiles)
    dirname = compound.inchikey
    if not os.path.exists(dirname):
        finder.target_smiles = compound.isomeric_smiles
        finder.tree_search()
        finder.build_routes()
        if finder.routes.images:
            os.mkdir(dirname)
            for n, image in enumerate(finder.routes.images):
                image.save(f"{dirname}/route{n:03d}.png")
            # Combine 4 images into a single image by ImageMagick
            images = glob.glob(f'{dirname}/route*.png')
            for i, imgs in enumerate(chunked(sorted(images), 3)):
                utils.concat_images(imgs, f"{dirname}/result{i}.png")
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
