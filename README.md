# discrete-hidden-markov-model
This is the implementation of Discrete Hidden Markov Model in C++. Mainly used to predict single-tandem-repetitive pattern in DNA after.

# Baum Welch EM Algorithm
The Baum-Welch Expectation-Maximization (oftenly called Re-Estimation) is a special type of unsupervised learning to cluster DNA. One example usage for G+C Rich pattern in DNA. The following implementation on C++ was derived from a paper by Rabiner in this [paper](https://www.ece.ucsb.edu/Faculty/Rabiner/ece259/Reprints/tutorial%20on%20hmm%20and%20applications.pdf), you may check that paper for the math background of the EM algorithm. All of the equations has been adapted to logarithmic operation in order to prevent the problem of **Vanishing Probability**.

# Viterbi Labeling (Clustering) Algorithm
The Viterbi Algorithm used to solve the secondary problem of Hidden Markov Model. This algorithm mainly used to reveal the *hidden* label for each corresponding nucleotides. In the area of G+C Rich Analysis, this labels represent the cluster in which the G and C nucleotides occurs.
