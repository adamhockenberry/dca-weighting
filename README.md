# dca-weighting
This repository is associated with the published manuscript titled: [Phylogenetic weighting does little to improve the accuracy of evolutionary coupling analyses](https://www.mdpi.com/1099-4300/21/10/1000). The code here should be sufficient to replicate all of the analyses presented in this manuscript.

The most important thing to note is that the `Data` and `Results` folders are not included here due to their large size/number of files. So in order to do more than look at / read through the code, users are required to [download the data archive from Zenodo](https://zenodo.org/record/3368652) (note that there are now two versions here, one from initial manuscript submission and one from the revised manuscript).

To construct all of the primary data, there are a few relevant steps and external programs that you'll have to call for each protein of interest. Here, I'll list the bash commands that I used for my system to do this.

First, I downloaded the PSICOV dataset from [here](http://bioinfadmin.cs.ucl.ac.uk/downloads/contact_pred_datasets/). (Note that this link may or may not always work and all of this data is available in archived form at the aforementioned Zenodo archive.)

You'll next want to go through the `clean_psicov_dataset.ipynb` notebook to do a bit of cleaning to the raw dataset and especially to create the random sample of sequences since I analyzed a maximum of 1001 sequences for each dataset. 

With the truncated/cleaned `.fasta` files, I used FastTree2 to create rough trees. I then ased IQTree to hopefully clean those rough trees up a bit in terms of branch lengths (and add that all of this is fully described in the [published](https://www.mdpi.com/1099-4300/21/10/1000) and [pre-print](https://www.biorxiv.org/content/10.1101/736173v1) versions of the manuscript as well as the markdown text in some of the notebooks).

Now after running these commands on all proteins and pushing some files around to make sure your directory structures match up you should be able to step through the `writing_weights.ipynb` and `comparing_weights.ipynb` files.

Since you've got the weights written to files now, it's time to run CCMPredPy which I did with some bash scripting. (note: I found that you'll want to process all the output matrices from this to deal with the meta-data contained in them).

And finally you should have all the files ready to go for the final analysis presented in `comparing_accuracies.ipynb`. 

