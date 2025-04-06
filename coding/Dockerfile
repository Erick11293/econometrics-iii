FROM jupyter/datascience-notebook:latest

USER ${NB_UID}

RUN pip install dbfread
RUN pip install geopandas
RUN pip install osmnx


USER root

RUN mkdir /home/compras
RUN cd /home/compras
WORKDIR /home/compras

RUN mkdir /home/compras/data
RUN mkdir /home/compras/src
RUN mkdir /home/compras/notebooks
RUN mkdir /home/compras/models


ENV NOTEBOOK_ARGS="--NotebookApp.token='' --NotebookApp.password=''"
COPY requirements.txt .



EXPOSE 8888