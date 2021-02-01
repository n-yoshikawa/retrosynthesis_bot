import argparse
import glob
import os

from more_itertools import chunked
from PIL import Image, ImageDraw
from aizynthfinder.aizynthfinder import AiZynthFinder

# https://stackoverflow.com/questions/30227466/
def concat_images(images, filename, padding=40):
    images = [Image.open(x) for x in images]
    widths, heights = zip(*(i.size for i in images))
    total_height = sum(heights) + (len(images)-1) * padding
    max_width = max(widths)
    new_img = Image.new('RGB', (max_width, total_height), (255, 255, 255))
    draw = ImageDraw.Draw(new_img)
    offset = 0
    for img in images:
        new_img.paste(img, (0, offset))
        offset += img.height + padding // 2
        draw.line(((0, offset), (new_img.width, offset)), fill=(0, 0, 0), width=2)
        offset += padding // 2
    new_img.save(filename)

def generate_images(smiles, out_dir, config):
    finder = AiZynthFinder(configfile=config)
    finder.stock.select("zinc")
    finder.expansion_policy.select("uspto")
    finder.filter_policy.select("uspto")

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    finder.target_smiles = smiles
    finder.tree_search()
    finder.build_routes()
    if finder.routes.images:
        for n, image in enumerate(finder.routes.images):
            image.save(f"{out_dir}/route{n:03d}.png")
        # Combine 4 images into a single image by ImageMagick
        images = glob.glob(f'{out_dir}/route*.png')
        for i, imgs in enumerate(chunked(sorted(images), 3)):
            concat_images(imgs, f"{out_dir}/result{i}.png")

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Random Retrosynthesis Bot')
    parser.add_argument('--smiles', required=True)
    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--config', default='/home/ubuntu/aizynthfinder/model/config.yml')
    args = parser.parse_args()

    generate_images(args.smiles, args.out_dir, args.config)
