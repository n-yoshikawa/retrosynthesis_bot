FROM continuumio/miniconda3
WORKDIR /root
RUN git clone https://github.com/MolecularAI/aizynthfinder.git
WORKDIR /root/aizynthfinder
RUN conda env create -f env-dev.yml
SHELL ["conda", "run", "-n", "aizynth-dev", "/bin/bash", "-c"]
RUN poetry install
RUN mkdir model && download_public_data model/
WORKDIR /root
RUN git clone https://github.com/n-yoshikawa/retrosynthesis_bot.git
WORKDIR /root/retrosynthesis_bot
RUN pip install tweepy pubchempy
