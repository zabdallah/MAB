# Multi- armed- Badit for Named Entity Recognition: 

This code runs multi armed bandit test (MAB) that have been used for Google analytics for evaluating the performance of Text mining tools (we focus on API based ones). The algorithm is described in this paper https://static.googleusercontent.com/media/research.google.com/en//pubs/archive/42550.pdf


MAB describes a hypothetical experiment where you face several slot machines ("one-armed bandits") with potentially different expected payouts. You want to find the slot machine with the best payout rate, but you also want to maximize your winnings. The algorithm maintains the posterior distribution of the payout rate of each arm, plays the arm in proportion to the probability that the given arm is optimal under the posterior, and then updates the posterior as new observations are made.

The algorithm stopping criteria is based on the probability that each variation beats the original. If we’re 70-80% sure that a variation beats the original then a winner has been found. 

The algorithm is briefly summarised as follows:
![alt text](https://github.com/zabdallah/MAB/blob/master/graph.png)

In bandit class, we first parse the command line for running commands that have a set of parameters such as methods to apply, what to extract and on which corpus. We assume the task here is NER (named entity extraction).  The tools are encapsulated in gapp files, that we have set up to call the APIs. The class uses some of GATE framework libraries (open source developed in university of Sheffield for text mining) https://gate.ac.uk/.

The maestro function in this class is **executeOnCorpus** . It first calls the **shuffleArray** method to randomise the corpus. All tools are applied on the first document (the seed). Then the performance measures are calculated for each. There are set of methods for calculating the statistics and performance metrics (**updateMetrics, calculateCopusStats, getMeasureValue, corpusQA**). Some of these have been overrode from GATE library  standard methods. The MAB bandit algorithm is implemented in the loop until convergence (which is checked with **converganceCheck** method returning Boolean **cont** ). The beta method generates beta sample based on the performance so far. Then, we choose the winner with the greatest value for the next iteration (argmax method). **Totable** and **WriteToFile** are for generating the output file in the write CSV format for visualisation.

The following is a sample graph based on the output file.
![alt text](https://github.com/zabdallah/MAB/blob/master/10sys.png)
