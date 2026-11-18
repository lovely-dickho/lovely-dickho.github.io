---
title: "[Econometrics] Maximum Likelihood Estimation"
date: 2026-11-18 16:24:00 +0900
categories: [Econometrics, Notes]
tags: [MLE, LR, Wald, LM, Fisher Information, CRLB]
math: true
toc: true
---

## MLE

우리는 이미 관측된 데이터(표본)를 가지고 있다. MLE는 “어떤 모수(Parameter)를 가정해야 지금 내 눈앞에 있는 이 데이터들이 관측될 확률이 가장 높아질까?”라는 질문에서 출발한다. 즉, 관측된 표본 데이터들이 나올 가능성(Likelihood)을 최대로 만들어주는 모수를 찾아내는 추정 방법이다.

#### 확률

모수(원인)가 주어졌을 때 특정 데이터(결과)가 나올 가능성이다. 예를 들어, 동전 앞면이 나올 확률이 0.5일 때 10번 던져서 앞면이 7번 나올 확률을 구하는 것이 여기에 해당한다.

#### 가능도

데이터(결과)가 이미 주어졌을 때, 이 데이터가 특정 모수(원인)로부터 나왔을 가능성이다. 예를 들어, 10번 던져서 앞면이 7번 나온 상황에서 동전 앞면이 나올 진짜 확률이 0.5일 가능성과 0.7일 가능성을 비교하는 것이 여기에 해당한다.

MLE는 바로 이 가능도를 극대화하는 수치를 최적의 추정량으로 선택하는 논리이다.

## MLE의 도출

확률변수 $$X$$가 모수 $$\theta$$를 가지는 확률밀도함수 $$f(x ; \theta)$$를 따른다고 가정한다. *iid*인 $$x_1, x_2, \cdots, x_n$$이 주어졌을 때, 이 표본들이 동시에 관측될 결합 확률밀도함수인 가능도 함수 $$L(\theta)$$는 개별 확률밀도함수의 곱으로 정의된다.

$$
L(\theta) = \prod_{i=1}^{n} f(x_i ; \theta)
$$

계산의 편의를 위해 가능도 함수 양변에 자연로그($$\ln$$)를 취해준다.

$$
\ln L(\theta) = \sum_{i=1}^{n} \ln f(x_i ; \theta)
$$

로그-가능도 함수를 미지의 모수 $$\theta$$에 대하여 편미분한 뒤, 그 기울기가 0이 되도록 하는 일계 조건 방정식을 세운다.

$$
\frac{\partial \ln L(\theta)}{\partial \theta} = 0
$$

### 정규분포에서의 최우추정량 도출 예시

확률변수 $$X$$가 평균이 $$\mu$$, 분산이 $$\sigma^2$$인 정규분포 $$N(\mu, \sigma^2)$$를 따른다고 가정한다. 확률밀도함수는 다음과 같다.

$$
f(x_i ; \mu, \sigma^2) = \frac{1}{\sqrt{2\pi\sigma^2}} \exp\left( -\frac{(x_i - \mu)^2}{2\sigma^2} \right)
$$

위 식을 전체 $$n$$개의 표본에 대해 곱하여 가능도 함수 $$L(\mu, \sigma^2)$$를 구성하고, 양변에 자연로그를 취해 $$\ln L(\mu, \sigma^2)$$를 구한다.

$$
\begin{aligned}
\ln L(\mu, \sigma^2) &= \sum_{i=1}^{n} \left[ -\frac{1}{2} \ln(2\pi) - \frac{1}{2} \ln(\sigma^2) - \frac{(x_i - \mu)^2}{2\sigma^2} \right] \\
&= -\frac{n}{2} \ln(2\pi) - \frac{n}{2} \ln(\sigma^2) - \frac{1}{2\sigma^2} \sum_{i=1}^{n} (x_i - \mu)^2
\end{aligned}
$$

위 식을 $$\mu$$에 대하여 편미분하고 0으로 놓는다.

$$
\frac{\partial \ln L}{\partial \mu} = - \frac{1}{2\sigma^2} \sum_{i=1}^{n} 2(x_i - \mu)(-1) = 0
$$

도출된 평균의 최우추정량은 표본평균과 일치한다.

$$
\hat{\mu}_{MLE} = \frac{1}{n} \sum_{i=1}^{n} x_i = \bar{x}
$$

### 분산(σ²)의 최우추정량

동일한 로그-가능도 함수를 이번에는 분산 $$\sigma^2$$ 자체를 하나의 변수로 보고 편미분하여 0으로 놓는다.

$$
\frac{\partial \ln L}{\partial (\sigma^2)} = -\frac{n}{2\sigma^2} + \frac{1}{2(\sigma^2)^2} \sum_{i=1}^{n} (x_i - \mu)^2 = 0
$$

양변에 $$2(\sigma^2)^2$$를 곱하여 식을 정리한다.

$$
-n\sigma^2 + \sum_{i=1}^{n} (x_i - \mu)^2 = 0 \quad \Longrightarrow \quad n\sigma^2 = \sum_{i=1}^{n} (x_i - \mu)^2
$$

위에서 구한 $$\hat{\mu}_{MLE}$$를 대입하여 $$\sigma^2$$에 대해 정리하면 분산의 최우추정량이 도출된다.

$$
\hat{\sigma}^2_{MLE} = \frac{1}{n} \sum_{i=1}^{n} (x_i - \bar{x})^2
$$

이는 OLS에서 구한 불편추정량 $$\hat{\sigma}^2 = \frac{1}{n-k}\sum e_i^2$$과 달리 분모가 $$n$$이므로 biased 추정량이다. ($$\mathbb{E}[\hat{\sigma}^2_{MLE}] = \frac{n-1}{n}\sigma^2 \neq \sigma^2$$)

## MLE in Matrix Form

#### 기본 가정 및 모형 설정

오차항에 대한 가정에 따라 $$u \sim N(0, \sigma^2 I_n)$$이 성립한다. 설명변수 행렬 $$X$$가 주어졌을 때, 종속변수 벡터 $$y$$의 조건부 분포는 $$y|X \sim N(X\beta, \sigma^2 I_n)$$이 된다.

#### Likelihood function

종속변수 $$y$$의 결합확률밀도함수 $$L(\beta, \sigma^2)$$는 다음과 같이 정의된다.

$$
L(\beta, \sigma^2) = (2\pi\sigma^2)^{-\frac{n}{2}} \exp\left( -\frac{1}{2\sigma^2} (y - X\beta)^\mathsf{T}(y - X\beta) \right)
$$

양변에 자연로그를 취하여 $$\ln L(\beta, \sigma^2)$$를 구한다.

$$
\ln L(\beta, \sigma^2) = -\frac{n}{2} \ln(2\pi) - \frac{n}{2} \ln(\sigma^2) - \frac{1}{2\sigma^2} (y - X\beta)^\mathsf{T}(y - X\beta)
$$

최우추정량 $$\hat{\beta}_{MLE}$$와 $$\hat{\sigma}^2_{MLE}$$를 구하기 위해 로그 우도함수를 각각의 파라미터로 편미분하여 0으로 둔다.

#### MLE for β

로그 우도함수에서 $$\beta$$와 관련된 항은 이차형식 $$(y - X\beta)^\mathsf{T}(y - X\beta)$$이다. 이를 전개하면 다음과 같다.

$$
(y - X\beta)^\mathsf{T}(y - X\beta) = y^\mathsf{T}y - 2\beta^\mathsf{T}X^\mathsf{T}y + \beta^\mathsf{T}X^\mathsf{T}X\beta
$$

로그 우도함수를 $$\beta$$로 편미분하면 다음과 같은 식을 얻는다.

$$
\frac{\partial \ln L}{\partial \beta} = -\frac{1}{2\sigma^2} \left( -2X^\mathsf{T}y + 2X^\mathsf{T}X\beta \right) = 0
$$

위 식을 정리하면 normal equation이 도출된다.

$$
X^\mathsf{T}X\beta = X^\mathsf{T}y \quad \Longrightarrow \quad \hat{\beta}_{MLE} = (X^\mathsf{T}X)^{-1}X^\mathsf{T}y
$$

&#8203;$$\beta$$의 최우추정량은 OLS 추정량과 동일하게 도출된다.

#### MLE for σ²

로그 우도함수를 분산 $$\sigma^2$$에 대해 편미분하면 다음과 같다.

$$
\frac{\partial \ln L}{\partial \sigma^2} = -\frac{n}{2\sigma^2} + \frac{1}{2(\sigma^2)^2} (y - X\beta)^\mathsf{T}(y - X\beta) = 0
$$

식을 정리하면 다음을 얻는다.

$$
n\sigma^2 = (y - X\beta)^\mathsf{T}(y - X\beta)
$$

여기서 $$\beta$$ 자리에 앞서 구한 최우추정량 $$\hat{\beta}_{MLE}$$를 대입하면, 잔차 벡터를 $$\hat{u} = y - X\hat{\beta}_{MLE}$$라 정의할 때 오차항 분산의 최우추정량은 다음과 같이 결정된다.

$$
\hat{\sigma}^2_{MLE} = \frac{(y - X\hat{\beta}_{MLE})^\mathsf{T}(y - X\hat{\beta}_{MLE})}{n} = \frac{\hat{u}^\mathsf{T}\hat{u}}{n}
$$

#### 결론

오차항이 정규분포를 따를 때, $$\beta$$의 최우추정량은 OLS 추정량과 완전히 일치한다. 반면, 분산 $$\sigma^2$$의 최우추정량은 분모가 자유도 수정을 거치지 않은 $$n$$이 되므로, OLS의 unbiased estimator (분모가 $$n-k$$) 와는 차이가 존재한다.

## MLE for the Uniform Distribution

연속균등분포 $$U(a, b)$$의 확률밀도함수는 다음과 같다. (단, $$a < b$$)

$$
f(x_i; a, b) = \frac{1}{b - a} \mathbf{1}_{\{a \le x_i \le b\}}
$$

우도함수 $$L(a, b)$$는 표본들의 지시함수 곱으로 전개된다.

$$
L(a, b | \mathbf{x}) = \frac{1}{(b - a)^n} \prod_{i=1}^n \mathbf{1}_{\{a \le x_i \le b\}}
$$

표본의 최소순서통계량을 $$x_{(1)} = \min(x_1, \dots, x_n)$$, 최대순서통계량을 $$x_{(n)}$$이라 할 때, 지시함수의 결합 조건은 $$a \le x_{(1)} \le x_{(n)} \le b$$이다.

$$
L(a, b | \mathbf{x}) =
\begin{cases}
\frac{1}{(b - a)^n}, & \text{if } a \le x_{(1)} \text{ and } b \ge x_{(n)} \\
0, & \text{otherwise}
\end{cases}
$$

#### 모수 (a, b)의 최적화

우도함수 $$L(a, b)$$를 최대화하는 것은 분모인 $$(b - a)^n$$을 최소화하는 것, 즉 구간의 길이 $$(b - a)$$를 최소화하는 것과 같다.

제약조건 $$a \le x_{(1)}$$ 및 $$b \ge x_{(n)}$$ 하에서 구간 $$(b - a)$$를 최소화하기 위한 필요충분조건은 다음과 같다.

**①** 하한 모수 $$a$$는 제약이 허용하는 가장 큰 값이어야 하므로, $$\hat{a} = x_{(1)}$$이다.

**②** 상한 모수 $$b$$는 제약이 허용하는 가장 작은 값이어야 하므로, $$\hat{b} = x_{(n)}$$이다.

따라서 모수 벡터 $$(a, b)$$의 최대우도추정량은 각각 표본의 최소 및 최대 순서통계량이다.

$$
\hat{a}_{MLE} = X_{(1)}, \quad \hat{b}_{MLE} = X_{(n)}
$$

#### 미분(FOC)을 통한 도출 불가능

미분을 본 모형에 적용할 수 없는 이유는 정칙조건(Regularity condition)이 위배되기 때문이다. 확률변수 $$X$$의 지지집합(Support)이 추정하고자 하는 모수 $$\theta$$에 의존하지 않아야 미분을 적용할 수 있다. 연속균등분포의 지지집합 $$[a, b]$$는 모수 자체에 의해 적분 구간의 경계가 직접 결정된다. 이로 인해 우도함수가 모수 공간의 경계점에서 불연속이 되어 미분이 불가능하므로, 통상적인 미분 연산자 대신 순서통계량을 활용한 경계점 최적화(Boundary optimization) 방식을 채택해야만 해를 도출할 수 있다.

## MLE for Dependent Data

표본 $$y_1, y_2, \cdots, y_T$$가 서로 독립적이지 않다면, 결합 확률밀도함수를 단순히 개별 확률밀도함수의 곱으로 나타낼 수 없다. 이 경우 확률의 연쇄 법칙(chain rule of probability)을 사용하여 조건부 확률의 곱으로 결합 확률밀도함수를 전개해야 한다.

$$
f(y_1, y_2, \cdots, y_T) = f(y_1) f(y_2 | y_1) f(y_3 | y_1, y_2) \cdots f(y_T | y_1, y_2, \cdots, y_{T-1})
$$

시계열 데이터에서 주로 나타나는 이와 같은 종속성은, 과거의 관측치가 주어졌을 때 현재 관측치가 나타날 조건부 확률을 극대화하는 방향으로 최우추정량을 도출하게 만든다.

가장 대표적인 종속 데이터 모형인 AR(1)을 생각해보자.

$$
y_t = \rho y_{t-1} + u_t
$$

여기서 오차항 $$u_t$$는 독립적이고 평균이 0, 분산이 $$\sigma^2$$인 정규분포 $$N(0, \sigma^2)$$를 따른다고 가정한다. 과거의 값 $$y_{t-1}$$이 주어졌을 때, 현재의 값 $$y_t$$는 평균이 $$\rho y_{t-1}$$이고 분산이 $$\sigma^2$$인 조건부 정규분포를 따른다.

$$
f(y_t | y_{t-1} ; \rho, \sigma^2) = \frac{1}{\sqrt{2\pi\sigma^2}} \exp\left( -\frac{(y_t - \rho y_{t-1})^2}{2\sigma^2} \right)
$$

Markov Property에 의해 $$y_t$$는 직전 과거인 $$y_{t-1}$$에만 직접적인 영향을 받으므로, 첫 번째 관측치 $$y_1$$이 주어졌을 때의 conditional likelihood function는 다음과 같이 표현된다.

$$
L(\rho, \sigma^2 | y_1) = \prod_{t=2}^{T} f(y_t | y_{t-1} ; \rho, \sigma^2)
$$

계산의 편의를 위해 양변에 자연로그를 취하여 conditional log-likelihood function을 도출한다.

$$
\begin{aligned}
\ln L(\rho, \sigma^2 | y_1) &= \sum_{t=2}^{T} \left[ -\frac{1}{2} \ln(2\pi) - \frac{1}{2} \ln(\sigma^2) - \frac{(y_t - \rho y_{t-1})^2}{2\sigma^2} \right] \\
&= -\frac{T-1}{2} \ln(2\pi) - \frac{T-1}{2} \ln(\sigma^2) - \frac{1}{2\sigma^2} \sum_{t=2}^{T} (y_t - \rho y_{t-1})^2
\end{aligned}
$$

위 식을 $$\rho$$에 대하여 편미분하고 0으로 놓아 1계 조건을 찾는다.

$$
\frac{\partial \ln L}{\partial \rho} = - \frac{1}{2\sigma^2} \sum_{t=2}^{T} 2(y_t - \rho y_{t-1})(-y_{t-1}) = 0 \quad \Longrightarrow \quad \frac{1}{\sigma^2} \sum_{t=2}^{T} (y_t y_{t-1} - \rho y_{t-1}^2) = 0
$$

이 방정식을 $$\rho$$에 대하여 정리하면 다음과 같이 최우추정량이 도출된다.

$$
\hat{\rho}_{MLE} = \frac{\sum_{t=2}^{T} y_t y_{t-1}}{\sum_{t=2}^{T} y_{t-1}^2}
$$

이는 조건부 최우추정법을 사용할 경우, 시계열 데이터에서도 OLS와 동일한 형태의 추정량이 도출됨을 보여준다.

## MLE와 검정의 직관

경제학자들이 많이 사용하는 test로는 LM test, LR test, White Test가 있다. 이 세 가지 test는 모두 그 아이디어가 비슷하다. 이들은 앞서 살펴본 MLE와 밀접한 관련이 있다. 우선 확률변수 $$W_i$$ 가 다음과 같다고 하자.

$$
W_1, W_2, W_3, \cdots, W_n \quad : \quad iid \text{ with density } f(W, \theta_0)
$$

우리의 관심사는 parameter $$\theta_0$$에 대한 제약이 성립하는가? 하는 것이다. 지금까지 우리는 t-test를 할 때에 우선 귀무가설을 $$H_0 : \theta_0 = \theta^*$$라고 두었다. 그리고 표본 평균 $$\overline{X}$$ 를 구해서 이것을 표준화하였다. 우리가 만약 귀무가설을 $$\theta_0 = 0$$ 이라고 두고 통계량을 구하면 다음과 같다.

$$
T = \frac{\overline{X} - 0}{\frac{s}{\sqrt{n}}} \sim t(n-1)
$$

이때 t-test는 이 값이 얼마나 0 에 가까운지를 보는 것이었다. 만약에 실제 parameter가 $$\theta_0 \neq 0$$ 이라면 t-value는 0에서 멀어지게 될 것이다. 다시 말해서 $$\theta_0$$ 의 일치추정량 $$\hat{\theta}$$은 $$\theta_0$$값에 converge 할 것이고, 만약 귀무가설이 맞는다면 $$\hat{\theta}$$과 $$\theta^*$$값 사이의 거리도 가까워질 것이다.

MLE를 통해 $$\hat{\theta}$$ 를 추정했을 때 둘 사이의 거리가 어떻게 측정되는지 보자. 우선 marginal PDF가 주어졌을 때, $$iid$$ 이면 joint PDF는 marginal PDF의 곱으로 표현할 수 있다. 계산의 편의를 위해 $$\ln$$을 취해 합으로 바꾸자.

$$
\text{Joint PDF} = L(\theta) = \prod_{i=1}^n f(W_i, \theta) \quad \longrightarrow \quad  \ln L(\theta) = \sum_{i=1}^n \ln f(W_i, \theta)
$$

이 joint PDF가 바로 likelihood function이다. $$\ln L(\theta)$$가 극대화되는 $$\theta$$를 찾는 것이 MLE 이다.

![MLE 그림](/assets/img/posts/mle-figure.png){: w="600" }

여기서 실제 $$\theta_0$$가 $$\theta^*$$라면 $$\theta^*$$와 $$\hat{\theta}$$ 사이의 거리도 가까울 것이고, 실제 $$\theta_0$$가 $$\theta^*$$가 아니라면 $$\theta^*$$와 $$\hat{\theta}$$ 사이의 거리는 멀어질 것이다. 즉, $$\hat{\theta}$$은 실제 $$\theta_0$$ 주위에 있기 때문에 $$\hat{\theta}$$와 $$\theta^*$$ 사이의 거리를 보고 test를 할 수 있는 것이다.

- **White Test**: 이 거리를 직접 측정
- **LR test**: $$\hat{\theta}$$와 $$\theta^*$$ 사이의 거리가 멀면 $$\ln L(\hat{\theta})$$과 $$\ln L(\theta^*)$$ 사이의 거리도 멀어진다는 점에 착안하여 함숫값의 거리를 측정
- **LM test**: 거리를 측정하는 대신 기울기 값을 보는 것이다. likelihood 값을 극대화 하는 $$\hat{\theta}$$ 이란 1계 조건에서 $$\hat{\theta}$$에서 likelihood function의 기울기가 0 이라는 것이다. 그러므로 $$\theta^*$$에서의 기울기 값이 얼마나 0 에 가까운지를 보면 $$\theta^*$$과 $$\hat{\theta}$$가 얼마나 가까운지를 알 수 있을 것이다.

LM test의 장점은 세 가지 test 중에서 가장 좋은 검정력을 가지고 있다는 것이다. 또한 White test나 LR test의 경우에는 가설검정을 위해서 비제약하의 값인 $$\hat{\theta}$$ 혹은 $$\ln L(\hat{\theta})$$과 제약하의 값 $$\theta^*$$ 혹은 $$\ln L(\theta^*)$$ 를 둘 다 구해야하지만, LM test 에서는 비제약하의 값은 항상 0 이므로 제약하의 값만 구하면 된다는 간편한 점이 있다.

## Likelihood Ratio Test (LR Test)

**무제약 모형과 제약 모형을 각각 추정하여, 두 모형이 산출하는 Max Likelihood의 높이 차이를 비교**한다.

**귀무가설과 통계량:**
무제약 추정량을 $$\hat{\theta}$$, 귀무가설 하의 제약 추정량을 $$\tilde{\theta}$$라고 할 때, 두 우도값의 비율에 자연로그를 취하고 $$-2$$를 곱하여 통계량을 구성한다.

$$
LR = -2 \ln \left( \frac{L(\tilde{\theta})}{L(\hat{\theta})} \right) = 2 \left( \ln L(\hat{\theta}) - \ln L(\tilde{\theta}) \right) \sim_{asy} \chi^2_q
$$

귀무가설 하에서 이 통계량은 점근적으로 제약식의 개수 $$q$$를 자유도로 하는 카이제곱 분포를 따른다.

**주로 쓰이는 상황:** 제약 모형과 무제약 모형 모두 수치적으로 추정하기 용이할 때 쓰이며, 모형 간의 적합도를 직접적으로 비교할 때 사용된다.

## Wald Test

오직 **무제약 모형 하나만 추정**하여, **무제약 추정량이 귀무가설이 설정한 값으로부터 수평적으로 얼마나 멀리 떨어져 있는가를 추정량의 분산으로 표준화하여 측정**한다.

### 스칼라 버전

하나의 특정 파라미터 $$\beta_k$$가 특정한 값 $$c$$와 같은지 검정한다 ($$H_0: \beta_k = c$$).
무제약 추정량 $$\hat{\beta}_k$$와 $$c$$의 거리를 분산으로 나눈다.

$$
W = \frac{(\hat{\beta}_k - c)^2}{\text{Var}(\hat{\beta}_k)} = t^2 \sim_{asy} \chi^2_1
$$

이는 통상적인 $$t$$ 통계량의 제곱과 같으며, 점근적으로 자유도가 1인 카이제곱 분포를 따른다.

### 행렬 버전

여러 개의 파라미터에 동시에 제약을 가하는 귀무가설 $$H_0: R\boldsymbol{\beta} = \mathbf{r}$$ 을 검정한다.
(여기서 $$R$$은 $$q \times k$$ 제약 행렬, $$\mathbf{r}$$은 $$q \times 1$$ 상수 벡터이다.)
무제약 추정량 벡터 $$\hat{\boldsymbol{\beta}}$$을 대입한 오차 벡터 $$(R\hat{\boldsymbol{\beta}} - \mathbf{r})$$와 이의 공분산 행렬을 이용하여 quadratic form을 구성한다.

$$
W = (R\hat{\boldsymbol{\beta}} - \mathbf{r})^\mathsf{T} [R \cdot \text{Var}(\hat{\boldsymbol{\beta}}) \cdot R^\mathsf{T}]^{-1} (R\hat{\boldsymbol{\beta}} - \mathbf{r}) \sim_{asy} \chi^2_q
$$

귀무가설 하에서 이 통계량은 점근적으로 제약식의 개수 $$q$$를 자유도로 하는 카이제곱 분포를 따른다.

**주로 쓰이는 상황:** 제약 모형을 추정하기 까다로울 때 유용하며, OLS의 개별 계수 $$t$$ 검정이나 다중 제약 $$F$$ 검정이 본질적으로 왈드 검정의 유한 표본 형태이다.

## Lagrange Multiplier Test (LM Test)

오직 **제약 모형 하나만 추정**하여, **제약 하에서 구한 추정량을 무제약 우도함수의 도함수(기울기)에 대입했을 때 그 기울기가 0에서 얼마나 떨어져 있는가**를 측정한다.

**귀무가설과 통계량:**
로그 우도함수의 스코어 함수를 $$S(\theta) = \partial \ln L(\theta) / \partial \theta$$, 피셔 정보 행렬을 $$I(\theta)$$라고 할 때, 제약 추정량 $$\tilde{\theta}$$를 평가한 기울기의 크기를 분산으로 표준화한다.

$$
LM = S(\tilde{\theta})^\mathsf{T} I(\tilde{\theta})^{-1} S(\tilde{\theta}) \sim_{asy} \chi^2_q
$$

**주로 쓰이는 상황:** 이미 추정된 기본 모형에 새로운 변수나 복잡한 동학을 추가했을 때 모형이 유의미하게 개선되는지 진단할 때 주로 쓰인다. 잔차 분석만으로 검정이 가능하여 다양한 진단 검정의 이론적 근간이 된다.

## Fisher Information Matrix and Cramér-Rao Lower Bound

**①** 스코어 벡터 $$S(\boldsymbol{\theta})$$의 기댓값은 $$\mathbf{0}$$이며, 그 공분산 행렬이 피셔 정보 행렬 $$\mathcal{I}(\boldsymbol{\theta})$$가 된다.

**②** 피셔 정보 행렬 $$\mathcal{I}(\boldsymbol{\theta})$$와 unbiased estimator $$\hat{\boldsymbol{\theta}}$$의 분산 하한은 서로 역행렬 관계에 있다. 즉, 로그 가능도 함수의 평균적인 곡률이 가파를수록 데이터가 참값 주변에 밀집해 있음을 의미하며, 이에 따라 추정 오차의 하한은 작아진다.

#### Mathematical Model

&#8203;$$k \times 1$$ 파라미터 벡터 $$\boldsymbol{\theta} = [\theta_1, \theta_2, \dots, \theta_k]^{\mathsf{T}}$$에 대하여, 로그 가능도 함수를 $$\ell(\boldsymbol{\theta}) \equiv \log L(\boldsymbol{\theta} \mid \mathbf{X})$$라 할 때 구체적인 행렬 형태는 다음과 같다.

Score Vector

$$
S(\boldsymbol{\theta}) = \nabla_{\boldsymbol{\theta}} \ell(\boldsymbol{\theta}) = \begin{bmatrix} \frac{\partial \ell}{\partial \theta_1} & \frac{\partial \ell}{\partial \theta_2} & \cdots & \frac{\partial \ell}{\partial \theta_k} \end{bmatrix}^{\mathsf{T}}
$$

Hessian Matrix

$$
\mathbf{H}_{\ell}(\boldsymbol{\theta}) = \begin{bmatrix}
\frac{\partial^2 \ell}{\partial \theta_1^2} & \frac{\partial^2 \ell}{\partial \theta_1 \partial \theta_2} & \cdots & \frac{\partial^2 \ell}{\partial \theta_1 \partial \theta_k} \\[8pt]
\frac{\partial^2 \ell}{\partial \theta_2 \partial \theta_1} & \frac{\partial^2 \ell}{\partial \theta_2^2} & \cdots & \frac{\partial^2 \ell}{\partial \theta_2 \partial \theta_k} \\[8pt]
\vdots & \vdots & \ddots & \vdots \\[8pt]
\frac{\partial^2 \ell}{\partial \theta_k \partial \theta_1} & \frac{\partial^2 \ell}{\partial \theta_k \partial \theta_2} & \cdots & \frac{\partial^2 \ell}{\partial \theta_k^2}
\end{bmatrix}
$$

Fisher Information

$$
\mathcal{I}(\boldsymbol{\theta}) = \text{Cov}(S(\boldsymbol{\theta})) = \mathbb{E}\left[ S(\boldsymbol{\theta}) S(\boldsymbol{\theta})^{\mathsf{T}} \right] = -\mathbb{E}\left[ \mathbf{H}_{\ell}(\boldsymbol{\theta}) \right]
$$

피셔 정보 행렬의 $$(i, j)$$번째 원소는 $$\mathcal{I}(\boldsymbol{\theta})_{i,j} = -\mathbb{E}\left[ \frac{\partial^2 \ell(\boldsymbol{\theta})}{\partial \theta_i \partial \theta_j} \right]$$로 정의된다.

#### 크라머-라오 하한(CRLB) 유도

임의의 불편추정량 $$\hat{\boldsymbol{\theta}}$$에 대하여 $$\mathbb{E}[\hat{\boldsymbol{\theta}}] = \boldsymbol{\theta}$$를 만족한다. 양변을 $$\boldsymbol{\theta}$$에 대해 편미분하면 적분과 미분의 교환 조건에 의해 다음의 공분산 관계가 유도된다.

$$
\text{Cov}(\hat{\boldsymbol{\theta}}, S(\boldsymbol{\theta})) = \mathbb{E}\left[ (\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}) S(\boldsymbol{\theta})^T \right] = \mathbf{I}_k
$$

일반화된 코시-슈바르츠 부등식(혹은 공분산 행렬의 사영 성질)에 의해, 임의의 두 확률벡터 $$X, Y$$에 대해 다음이 성립한다.

$$
\text{Cov}(X) \succeq \text{Cov}(X,Y)\text{Cov}(Y)^{-1}\text{Cov}(Y,X)
$$

여기에 $$X = \hat{\boldsymbol{\theta}}$$, $$Y = S(\boldsymbol{\theta})$$를 대입하면 최종적인 크라머-라오 부등식을 얻는다.

$$
\text{Cov}(\hat{\boldsymbol{\theta}}) \succeq \mathbf{I}_k \cdot \mathcal{I}(\boldsymbol{\theta})^{-1} \cdot \mathbf{I}_k = \mathcal{I}(\boldsymbol{\theta})^{-1}
$$
