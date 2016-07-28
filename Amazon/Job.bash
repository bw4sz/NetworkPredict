#!/bin/bash
#cd into occupy, if directory doesn't exist, kill the run.
git git@github.com:bw4sz/NetworkPredict.git --depth 1

cd NetworkPredict||sudo halt

#git pull to make sure we are at HEAD
git pull

#make new branch
#name it the instance ID
iid=$(ec2metadata --instance-id)

git checkout -b $iid

#render script
Rscript -e "rmarkdown::render('Observed2m.Rmd')" &> run.txt
Rscript -e "rmarkdown::render('Observed2m_model.Rmd')" &> run.txt

#push results
git add --all
git commit -m "ec2 run complete"
git push -u origin $iid

#kill instance
sudo halt
