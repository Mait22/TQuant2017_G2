### <span style="color: #84BD93">What we did </span>
Anlogous to the paper, we first set initial values for the parameter of the
probability functions: 
$$
M1: \pi_{i} = (1 + t_{i})^{-a}\\
M2: \pi_i = (b + t_{i})^{-a}\\
M3: \pi_i = (1 + b*t_{i})^{-a}
$$
Then data ist generated using a binomial distribution with $n = 50$, and using
the pi from the respective model, as well as $t = (.1, 2.1, 4.1, 6.1, 8.1)$.
With this data, maxmimum likelihood estimations were produced using the
likelihood $$L: \binom{n}{y_{i}} * \pi_{i}^{y_{i}} * (1 - \pi_{i})^{n - y_{i}}. $$

These estimation then were used as a basis for addition data generation. A
sampling error was introduced. In this case we used a proportion of the real
probability as well as a complementary proportion of a uniform distribution of
 values between 0 and 1. 

Having thus generated data from each model, each model then was fitted to this
generated data. For each fit GoF indices were computed (RMSE, PVAF) as well
as the AIC $(-2(loglike) + 2k)$ and BIC $(-2(loglike) +
k(ln(n)))$. This step was repeated N times, and the proportion of
choosing the respective model given a indice was then computed for each
dataset. 

Given the results from the paper, RMSE and PVAf should prove an ineffective
way in which the proper generator model will be selected, because they will
always go for the model with the best fit towards the data. Using AIC, or BIC
for example, should allow the user to not only fit the observed data, but also
select the model which might predict future data better. 

### <span style="color: #84BD93"> Glossary, from Myung, 2001</span>
  + Complexity: the property of a model that enables it to fit diverse patterns of data; it is the
flexibility of a model. Although the number of parameters in a model and its functional form can
be useful for gauging its complexity, a more accurate and intuitive measure is the number of
distinct probability distributions that the model can generate by varying its parameters over
their entire range. 

  + Functional form: the way in which the parameters $(\beta)$ and data (x) are combined in a model’s
equation: $y = \beta x$ and $y = \beta + x$ have the same number of parameters but different functional forms
(multiplicative versus additive). 

  + Generalizability: the ability of a model to fit all data samples generated by the same cognitive
process, not just the currently observed sample (i.e. the model’s expected GOF with respect to
new data samples). Generalizability is estimated by combining a model’s GOF with a measure of
its complexity. 

  + Goodness of fit (GOF): the precision with which a model fits a particular sample of observed
data. The predictions of the model are compared with the observed data. The discrepancy
between the two is measured in a number of ways, such as calculating the root mean squared
error between them. 

  + Overfitting: the case where, in addition to fitting the main trends in the data, a model also fits the
microvariation from this main trend at each data point. 

  + Parameters: variables in a model’s equation that represent mental constructs or processes;
they are adjusted to improve a model’s fit to data. For example, in the model $y = \beta x$, $\beta$ is a
parameter. 