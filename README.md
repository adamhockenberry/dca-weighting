# dca-weighting
This repository is for a forthcoming unsubmitted manuscript. This README will be updated continuously until publication. 

The most important thing to note is that the `Data` and `Results` folders are not included here due to size restrictions. So in order to do more than look at / read through the code, users are required to [download the data archive from Zenodo(but this does not yet work)]().

To construct all of the primary data, there are a few relevant steps and external programs that you'll have to call for each protein of interest. Here, I'll list the bash commands that I used for my system to do this.

First, I downloaded the PSICOV dataset from [here](http://bioinfadmin.cs.ucl.ac.uk/downloads/contact_pred_datasets/). (Note that this link may or may not always work and all of this data is available in archived form at the Zenodo archive.)

You'll next want to go through the `clean_psicov_dataset.ipynb` notebook to do a bit of cleaning to the raw dataset and especially to create the random sample of sequences since I analyzed a maximum of 1001 sequences for each dataset. 

With the truncated/cleaned `.fasta` files, I used FastTree2 to create rough trees as follows:

I then used IQTree to hopefully clean those rough trees up a bit in terms of branch lengths:

Now after running these commands on all proteins and pushing some files around to make sure your directory structures match up you should be able to step through the `writing_weights.ipynb` and `comparing_weights.ipynb` files.

Since you've got the weights written to files now, it's time to run CCMPredPy which I did as follows:

You'll want to process all the output matrices from this to deal with the meta-data contained in them:

And finally you should have all the files ready to go for the final analysis presented in `comparing_accuracies.ipynb`. 

