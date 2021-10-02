# DWOP

This is a MATLAB implementation of the article "Online hotel rating prediction through a dynamic weighted ordered probit model" written by Myrsini Ntemi and Constantine Kotropoulos and published in Digital Signal Processing, Elsevier, 2021.


# Abstract
Hotel reputation might volatile through time due to improvement or deterioration of accommodation facilities, management strategies changes, etc. A recommender system should take this information into consideration in order to provide  updated recommendations. The vast majority of recommender systems rely on matrix factorization techniques without explicitly model hotel reputation evolution through time. In real life, recommender systems face another challenge, namely data sparsity, which affects strongly prediction accuracy. In this paper, a dynamic weighted ordered probit model is integrated with collaborative Kalman filtering for online hotel rating prediction.   The model learns  a dynamically evolving optimal standard deviation via a Newton-Raphson type update, modeling the hotel reputation evolution. An accept-reject sampling scheme is applied to compute on the fly efficiently the expected value of an auxiliary truncated normal random variable utilized to predict hotel ratings based  on past ratings, which fill in a short window.  Such a window can always be filled in, confronting thus data sparsity and scalability issues. Experiments have demonstrated that the proposed model outperforms state-of-the-art techniques. Performance gains are attested to be statistically significant.


If you find this project useful, please cite:

@article{ntemi2020dynamic,
  title={Online hotel rating prediction through a dynamic weighted ordered probit model},
  author={Ntemi, Myrsini and Kotropoulos, Constantine},
  journal={Digital Signal Processing},
  volume={xxx},
  pages={xxxx},
  year={2021xxx?},
  publisher={Elsevier}
}


