import argparse
import sys

from DECIMER import predict_SMILES

parser = argparse.ArgumentParser('Retrosynthesis Bot')
parser.add_argument('--image')
args = parser.parse_args()

try:
    smiles = predict_SMILES(args.image)
except Exception as e:
    print(e, file=sys.stderr)
else:
    print(smiles)
