import os
import sys
import glob
import traceback
import subprocess
from more_itertools import chunked

import tweepy
from aizynthfinder.aizynthfinder import AiZynthFinder
import pubchempy as pcp
import random

if __name__ == "__main__":
    # Set Twitter API Access Tokens
    consumer_key = "XXXXXXXXXXXXXXXXXXXXXXXXX"
    consumer_secret = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    access_token = "XXXXXXXXXXXXXXXXXXX-XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    access_secret = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"

    auth = tweepy.OAuthHandler(consumer_key, consumer_secret)
    auth.set_access_token(access_token, access_secret)
    api = tweepy.API(auth)

    finder = AiZynthFinder(configfile="/home/ubuntu/network/config.yml")
    finder.stock.select("zinc")
    finder.expansion_policy.select("uspto")
    finder.filter_policy.select("uspto")

    while True:
        cid = random.randint(10000, 1000000)
        compound = pcp.get_compounds(cid, 'cid')[0]
        print(compound.isomeric_smiles)
        print(compound.synonyms)
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
                    arg = ' '.join(imgs)
                    proc = subprocess.run(f"convert -append -gravity south -splice 0x40 {arg} {dirname}/result{i}.png > /dev/null", shell=True, text=True)
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
