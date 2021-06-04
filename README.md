# uniprod
Why sample mean of products of n Uniform(0,1)Uniform(0,1) distributed independant random numbers is neither 1/2^{n}1/2 
n
  nor 1/e^{n}1/e 
n
 .
I happened to notice that when multiplying independant Uniform(0,1) random variables, the product goes to zero after around 745 variables. Smallest non zero 64-bit floating point number is 2^-1074 and one could think that going to zero would need a product of around 1074 variables. Assuming that E[x_1 \cdot x_2 \cdot...\cdot x_n] = 1/2^{n}E[x 
1
​	
 ⋅x 
2
​	
 ⋅...⋅x 
n
​	
 ]=1/2 
n
 . Also sampling enough of products of length n, the sample mean is < 1/2^{n}<1/2 
n
 , if n > about 20.

Below I have tried to deduce E[x_1 \cdot x_2 \cdot...\cdot x_n]E[x 
1
​	
 ⋅x 
2
​	
 ⋅...⋅x 
n
​	
 ] from distributions connected to Uniform(0, 1) variables. But the problem still remains: why sampled means do not fit to expected values from distributions.

Expectation of independant variables
When two random variables are statistically independent, the expectation of their product is the product of their expectations (e.g. Wikipedia).

x_i \sim U(0, 1)
x 
i
​	
 ∼U(0,1)
E[x_i] = 1/2
E[x 
i
​	
 ]=1/2
x = \prod_{i=1}^{n}x_i
x= 
i=1
∏
n
​	
 x 
i
​	
 
From the independance assumption we get the first nominee for E[x]

E[x] = \prod_{i=1}^{n} E[x_i] = 1/2^{n}
E[x]= 
i=1
∏
n
​	
 E[x 
i
​	
 ]=1/2 
n
 
Expected value for geometric mean of x (StackExchange) is

E[x^{1/n}] = (1 + 1/n)^{-n}
E[x 
1/n
 ]=(1+1/n) 
−n
 
n \to \infty \Rightarrow (1 + 1/n)^{-n} \rightarrow 1/e
n→∞⇒(1+1/n) 
−n
 →1/e
n \gtrapprox 100 \Rightarrow E[x^{1/n}] \cong 1/e
n⪆100⇒E[x 
1/n
 ]≅1/e
Random sampled geometric means go nicely near to (1 + 1/n)^{-n}(1+1/n) 
−n
 . If E[x^{1/n}] \cong 1/eE[x 
1/n
 ]≅1/e, then E[x]E[x] generally cannot be much bigger than 1/e^{n}1/e 
n
 ?
Uniform Product Distribution
(WolframAlpha, Wiki)

\gamma (n, k)γ(n,k) is the lower incomplete gamma function and \gamma (n, k) / (n-1)!γ(n,k)/(n−1)! is Gamma(n, k) distribution function. All integrals are from WolframAlpha.

f_u(x;n) = (-log(x))^{n-1} / (n-1)!
f 
u
​	
 (x;n)=(−log(x)) 
n−1
 /(n−1)!
E[x] = \int_{0}^{1} x\cdot f_u(x;n) dx
E[x]=∫ 
0
1
​	
 x⋅f 
u
​	
 (x;n)dx
E[x] = 1/2^{n}\cdot \gamma (n, 0) / (n-1)!
E[x]=1/2 
n
 ⋅γ(n,0)/(n−1)!
\gamma (n, 0) / (n-1)! = 1 \Rightarrow E[x]= 1/2^{n}
γ(n,0)/(n−1)!=1⇒E[x]=1/2 
n
 
Gamma and Exponential distributions
y = log(x) = \sum_{i=1}^{n} log(x_i)
y=log(x)= 
i=1
∑
n
​	
 log(x 
i
​	
 )
-log(x_i) \sim Exp(1) \Rightarrow E[y] = \sum_{i=1}^{n} E[log(x_i)] = -n
−log(x 
i
​	
 )∼Exp(1)⇒E[y]= 
i=1
∑
n
​	
 E[log(x 
i
​	
 )]=−n
-y \sim Gamma(n, 1) \Rightarrow E[y] = -n
−y∼Gamma(n,1)⇒E[y]=−n
E[x] = e^{E[y]} = 1/e^{n}?
E[x]=e 
E[y]
 =1/e 
n
 ?
Poisson distribution
Expected number of variables x_ix 
i
​	
  needed in the product to reach value e^{-\lambda}e 
−λ
  is distributed Poisson(\lambda)Poisson(λ). D. Knuth has even given an algorithm to generate Poisson distributed random numbers by multiplying U(0, 1) random numbers.

n \sim Poisson(\lambda)
n∼Poisson(λ)
f_p(n; \lambda) = \lambda^{n-1}\cdot e^{-\lambda} / (n-1)!
f 
p
​	
 (n;λ)=λ 
n−1
 ⋅e 
−λ
 /(n−1)!
The Poisson mass probability function can also be expressed by Gamma density function, where the normal parameters have changed place and role.

R: dgamma(shape=n, x=\lambda) = dpois(n-1, \lambda)
R:dgamma(shape=n,x=λ)=dpois(n−1,λ)
Mean number of x_ix 
i
​	
 's to reach a limit z is -log(z). Reaching the smallest 64-bit FP number 2^{-1074}2 
−1074
 , -log(2^{-1074}) = log(2) \cdot 1074 \cong 744.44−log(2 
−1074
 )=log(2)⋅1074≅744.44 variables are needed. The probability that a product of 1074 variables does't go to zero is 4e-30.

For fixed number n we can integrate the product mean by integrating mean of limit values e^{-\lambda}e 
−λ
  over \lambdaλ range 0 - \infty0−∞.

E[x] = \int_{0}^{\infty}e^{-\lambda}\cdot f_p(n,\lambda) d\lambda = 1/2^{n}\cdot \gamma (n, 0) / (n-1)!
E[x]=∫ 
0
∞
​	
 e 
−λ
 ⋅f 
p
​	
 (n,λ)dλ=1/2 
n
 ⋅γ(n,0)/(n−1)!
\gamma (n, 0) / (n-1)! = 1 \Rightarrow E[x] = 1/2^{n}
γ(n,0)/(n−1)!=1⇒E[x]=1/2 
n
 
Uniform Product Distribution and Poisson distribution give the same integral and result 1/2^{n}1/2 
n
 .

Here we have a kind of controversy: general mean of a products of length n is 1/2^{n}1/2 
n
  and a product of lenth n reaching limit e^{-n}e 
−n
  has mean 1/e^{n}1/e 
n
  and n is the expected number of variables for reaching the limit.

This is off topic, but here is the nice Knuth's Poisson random number generator in R.

# Knuth's algorithm in R
rpoisson <- function(lambda) {
    L <- exp(-lambda)
    k <- 0
    p <- 1
    while (TRUE) {
        k <- k + 1
        p <- p * runif(1)
        if (p < L) {
            return(k - 1)
        }
    }
}
Results of simulations
It is very clear that E[x_1 \cdot x_2] = 0.25E[x 
1
​	
 ⋅x 
2
​	
 ]=0.25 and E[x_1 \cdot x_2\cdot x_3] = 0.125E[x 
1
​	
 ⋅x 
2
​	
 ⋅x 
3
​	
 ]=0.125. Or at least differences are hard to detect in presence of numerical noise from FP arithmetic and possible very small biases in pseudo random number generators.

The sampling was mostly made by R with 64-bit FP arithmetic and checked also by C/C++ code with 128-bit arithmetic. 64-bits is enough for n \lessapprox 700n⪅700 and 128-bits for n \lessapprox 11000n⪅11000. For running \bar{X_n } 
X 
n
​	
 
ˉ
​	
  to the limit 1/e^{n}1/e 
n
  a special "small number multiplicator" R-function was coded and used. Also more than one uniform random number generator was used.

Let \bar{X_n } = 
X 
n
​	
 
ˉ
​	
 = arithmetic mean of sampled products of lenght n.

For n > 50 it comes immediately clear that 1/2^{n} < \bar{X_n } < 1/e^{n}1/2 
n
 < 
X 
n
​	
 
ˉ
​	
 <1/e 
n
 .

For n \lessapprox 100n⪅100 \bar{X_n} 
X 
n
​	
 
ˉ
​	
  keeps near 1/2^{n}1/2 
n
  and when n increases \bar{X_n} 
X 
n
​	
 
ˉ
​	
  very slowly converges towards 1/e^{n}1/e 
n
 . Actually reaching 1/e^{n}1/e 
n
  needs n > 100000.

For n = 600 \bar{X_n} 
X 
n
​	
 
ˉ
​	
  is quite middle of 1/2^{n}1/2 
n
  and 1/e^{n}1/e 
n
 .

Sampled geometric means and log(product) sets to their distribution values (1 + 1/n)^{-n}(1+1/n) 
−n
  and -n.

Question: What is happening and why?
