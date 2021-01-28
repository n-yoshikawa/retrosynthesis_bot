import os
import glob
import traceback

from PIL import Image, ImageDraw
from more_itertools import chunked

import tweepy
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from aizynthfinder.aizynthfinder import AiZynthFinder

# https://stackoverflow.com/questions/30227466/
def concat_images(images, filename):
    images = [Image.open(x) for x in images]
    widths, heights = zip(*(i.size for i in images))
    total_height = sum(heights)
    max_width = max(widths)
    new_img = Image.new('RGB', (max_width, total_height), (255, 255, 255))
    draw = ImageDraw.Draw(new_img)
    offset = 0
    for img in images:
        new_img.paste(img, (0, offset))
        offset += img.height + 20
        draw.line(((0, offset), (new_img.width, offset)), fill=(0, 0, 0), width=2)
        offset += 20
    new_img.save(filename)

finder = AiZynthFinder(configfile="/home/ubuntu/network/config.yml")
finder.stock.select("zinc")
finder.expansion_policy.select("uspto")
finder.filter_policy.select("uspto")
    
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
        dirname = Chem.inchi.MolToInchiKey(m)
        # if dir does not exist (i.e. first run), run retrosynthesis
        if not os.path.exists(dirname):
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
                finder.target_smiles = smiles
                finder.tree_search()
                finder.build_routes()
                for n, image in enumerate(finder.routes.images):
                    image.save(f"{dirname}/route{n:03d}.png")
                # Concat 3 images into a single image by Pillow
                images = glob.glob(f'{dirname}/route*.png')
                for i, imgs in enumerate(chunked(sorted(images), 3)):
                    concat_images(imgs, f"{dirname}/result{i}.png")
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
