#!/bin/bash

echo "facebook"
./ldf -i ../data/facebook.adj -o ../data/facebookComm.txt
echo "wiki"
./ldf -i ../data/wiki.adj -o ../data/wikiComm.txt
echo "epinions"
./ldf -i ../data/epinions.adj -o ../data/epinionsComm.txt
echo "dblp"
./ldf -i ../data/dblp.adj -o ../data/dblpComm.txt
echo "pokec"
./ldf -i ../data/pokec.adj -o ../data/pokecComm.txt
