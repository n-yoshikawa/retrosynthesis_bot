from PIL import Image, ImageDraw

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
