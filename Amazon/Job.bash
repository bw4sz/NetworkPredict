#!/bin/bash
#cd into repo, if directory doesn't exist, kill the run.
git clone git@github.com:bw4sz/NetworkPredict.git --depth 1

cd NetworkPredict||sudo halt

#git pull to make sure we are at HEAD
git pull

#make new branch
#name it the instance ID
iid=$(ec2metadata --instance-id)

git checkout -b $iid

#render script
Rscript -e "rmarkdown::render('Observed2m.Rmd')"
Rscript -e "rmarkdown::render('Observed2m_model.Rmd')" &> run.txt

#push results
git add --all
git commit -m "ec2 run complete"
git push -u origin $iid

#kill instance
sudo halt
