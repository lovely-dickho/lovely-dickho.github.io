---
title: "[Econometrics] Instrumental Variable"
date: 2026-11-18 16:23:00 +0900
categories: [Econometrics, Notes]
tags: [IV, 내생성, 도구변수, 약한 도구변수]
math: true
toc: true
---

## Introduction

단순회귀모형을 가정하자. 여기서 설명변수 $$x_i$$가 오차항 $$u_i$$와 상관관계를 가져 Endogeneity가 발생했다고 가정한다.

$$
y_i = \beta_0 + \beta_1 x_i + u_i \quad Cov(x_i, u_i) \neq 0
$$

이때, 다음 두 조건을 만족하는 도구변수 $$z_i$$를 도입한다.

**① Relevance:** $$Cov(z_i, x_i) \neq 0$$

**② Exogeneity:** $$Cov(z_i, u_i) = 0 \implies \mathbb{E}[z_i u_i] = 0$$

외생성 가정 $$\mathbb{E}[z_i u_i] = 0$$에 $$u_i = y_i - \beta_0 - \beta_1 x_i$$를 대입하여 표본 모멘트 조건을 구하고, 편차 형태로 나타내면 스칼라 형태의 IV 추정량이 도출된다.

$$
\hat{\beta}_{1,IV} = \frac{\sum_{i=1}^{n} (z_i - \bar{z})(y_i - \bar{y})}{\sum_{i=1}^{n} (z_i - \bar{z})(x_i - \bar{x})} = \frac{Cov(z,y)}{Cov(z,x)}
$$

혹은 다음과 같은 방법으로도 IV 추정량을 구할 수 있다. 회귀모형식에, 도구변수 $$z_i$$와 공분산을 취하자.

$$
\begin{aligned}
Cov(y_i, z_i) &= Cov(\beta_0 + \beta_1 x_i + u_i , z_i) \\
&= Cov(\beta_0, z_i) + Cov(\beta_1 x_i, z_i) + Cov(u_i, z_i) \quad (\because \text{Linearity of Covariance}) \\
&= \beta_1 Cov( x_i, z_i) \\
\therefore \beta_1 &= \frac{Cov(z,y)}{Cov(z,x)}
\end{aligned}
$$

공분산의 선형성, 도구변수의 외생성 가정을 통해 식을 정리하면 $$\beta_1$$ 얻을 수 있다.

## 불편성, 일치성, 효율성

IV 추정량은 소표본에서 bias가 발생한다.

$$
\mathbb{E}[\hat{\beta}_{1,IV}] = \beta_1 + \mathbb{E}\left[ \frac{\sum_{i=1}^{n} (z_i - \bar{z}) u_i}{\sum_{i=1}^{n} (z_i - \bar{z})(x_i - \bar{x})} \right] \neq \beta_1
$$

표본의 크기 $$n \to \infty$$일 때 IV 추정량은 모수에 확률적으로 수렴한다.

$$
\text{plim} (\hat{\beta}_{1,IV}) = \beta_1 + \frac{\text{plim} \left( \frac{1}{n} \sum_{i=1}^{n} (z_i - \bar{z}) u_i \right)}{\text{plim} \left( \frac{1}{n} \sum_{i=1}^{n} (z_i - \bar{z})(x_i - \bar{x}) \right)}
$$

대수의 법칙(LLN)과 Slutsky's Thm에 의해 극한은 모집단의 공분산으로 수렴한다.

$$
\text{plim} (\hat{\beta}_{1,IV}) = \beta_1 + \frac{Cov(z, u)}{Cov(z, x)}
$$

도구변수의 외생성 가정에 의해 $$Cov(z, u) = 0$$이므로,

$$
\text{plim} (\hat{\beta}_{1,IV}) = \beta_1
$$

따라서 IV 추정량은 **일치 추정량**이다.

## Asymptotic Variance

오차항의 동분산성 및 자기상관이 없음을 가정하자.

$$
Var(u_i) = \sigma^2, \quad Cov(u_i, u_j) = 0 \quad (i \neq j)
$$

도구변수 추정량은 편차 형태로 다음과 같이 얻었다.

$$
\hat{\beta}_{IV} = \beta_1 + \frac{\sum (z_i - \bar{z}) u_i}{\sum (z_i - \bar{z})(x_i - \bar{x})}
$$

설명변수와 도구변수를 고정된 상수로 취급하여 양변에 분산을 취하고 정리하면, 다음을 얻는다.

$$
\begin{aligned}
Var(\hat{\beta}_{IV}) &= \frac{\sum (z_i - \bar{z})^2 Var(u_i)}{\left[ \sum (z_i - \bar{z})(x_i - \bar{x}) \right]^2} \\
&= \frac{\sigma^2 \sum (z_i - \bar{z})^2}{\left[ \sum (z_i - \bar{z})(x_i - \bar{x}) \right]^2} \quad (\because \text{Homoskedasticity Assumption})
\end{aligned}
$$

위 식의 분자와 분모를 각각 $$\sum (z_i - \bar{z})^2 \sum (x_i - \bar{x})^2$$ 로 나누자.

$$
Var(\hat{\beta}_{IV}) = \frac{ \frac{\sigma^2 \sum (z_i - \bar{z})^2}{\sum (z_i - \bar{z})^2 \sum (x_i - \bar{x})^2} }{ \frac{\left[ \sum (z_i - \bar{z})(x_i - \bar{x}) \right]^2}{\sum (z_i - \bar{z})^2 \sum (x_i - \bar{x})^2} }
$$

분자에서는 $$\sum (z_i - \bar{z})^2$$ 이 약분되어 OLS의 분산 형태가 도출된다. 분모는 도구변수 $$z$$와 설명변수 $$x$$ 간의 표본 상관계수 제곱($$r_{xz}^2$$)이 된다.

$$
\begin{aligned}
\text{분자} &= \frac{\sigma^2}{\sum (x_i - \bar{x})^2} = Var(\hat{\beta}_{OLS}) \\
\text{분모} &= r_{xz}^2
\end{aligned}
$$

식을 정리하면 최종적으로 다음을 얻는다.

$$
Var(\hat{\beta}_{IV}) = \frac{Var(\hat{\beta}_{OLS})}{r_{xz}^2}
$$

동분산성 가정 하에서 IV 추정량의 점근적 분산은 OLS보다 크다.

$$
Var(\hat{\beta}_{1,IV}) \approx \frac{\sigma^2}{SST_x \cdot R_{x,z}^2} = \frac{Var(\hat{\beta}_{OLS})}{r_{xz}^2}
$$

일반적으로 상관계수의 제곱 $$r_{xz}^2 \le 1$$이므로, **IV 추정량의 분산은 항상 OLS 추정량의 분산보다 크다** (효율성이 떨어진다). 도구변수 $$z$$가 $$x$$를 약하게 설명할수록($$r_{xz}^2$$가 작을수록) 분산은 폭발적으로 증가하며, 이를 ‘약한 도구변수 문제’라고 한다.

## IV와 Proxy의 차이

연구에 필요한 핵심 변수(예: 능력, $$x^*$$)를 데이터로 관측할 수 없는 omitted variable 상황일 때, 이를 대신하는 관측 가능한 변수($$x$$)를 proxy variable라고 한다. 관측 불가능한 변수 $$x^*$$가 포함된 진짜 모델이 다음과 같다고 하자.

$$
Y = \beta_0 + \beta_1 x^* + \gamma W + u
$$

&#8203;$$x^*$$와 우리가 관측할 수 있는 proxy $$x$$의 관계 설정:

$$
x^* = \delta_0 + \delta_1 x + v
$$

이를 원래 식에 대입하여 실제로 추정하게 되는 방정식:

$$
Y = (\beta_0 + \beta_1 \delta_0) + (\beta_1 \delta_1) x + \gamma W + (\beta_1 v + u)
$$

proxy $$x$$가 $$x^*$$를 잘 설명하려면, 오차항 $$v$$가 다른 통제변수들($$W$$)이나 원래 방정식의 오차항 $$u$$와 상관관계가 없어야 한다. 즉, proxy는 원래의 오차항과 독립적이어야 편의를 최소화할 수 있다.

도구변수와 proxy의 차이를 정리하면 다음과 같다.

#### IV

변수 $$X$$를 **관측할 수는 있지만**, 오차항과 상관관계(내생성)가 있을 때 이를 해결하기 위해 외부에서 투입하는 ‘외생적 충격’ 변수 $$Z$$.

#### Proxy

변수 $$X^*$$를 **아예 관측할 수 없을 때**, 모형을 추정하기 위해 꿩 대신 닭으로 투입하는 ‘대체’ 변수 $$X$$.

## 도구변수의 강도 체크

도구변수 $$Z$$는 내생변수 $$X$$와 충분히 높은 상관관계를 가져야 한다($$Cov(Z, X) \neq 0$$). 이 상관관계가 약하면 ‘약한 도구변수’ 문제가 발생하여, 추정량의 분산이 커지고 finite sample bias가 심각해진다.

#### 도구변수가 1개일 때

내생변수 $$X_t$$를 도구변수 $$Z_t$$에 대해 회귀하는 1단계 회귀를 설정한다.

$$
X_t = \pi_0 + \pi_1 Z_t + v_t
$$

여기서 귀무가설 $$H_0: \pi_1 = 0$$에 대한 $$F$$-통계량(또는 $$t$$-통계량의 제곱)을 확인한다. 일반적으로 1단계 $$F$$-통계량이 10보다 작으면 도구변수가 약하다고 판단한다.

#### 도구변수가 여러 개일 때

&#8203;$$k$$개의 도구변수 $$Z_{1t}, \dots, Z_{kt}$$가 있을 때의 1단계 회귀식은 다음과 같다.

$$
X_t = \pi_0 + \pi_1 Z_{1t} + \dots + \pi_k Z_{kt} + v_t
$$

이때 귀무가설 $$H_0: \pi_1 = \dots = \pi_k = 0$$에 대한 결합 유의성을 $$F$$-검정으로 확인한다. 마찬가지로 $$F$$-통계량이 10을 넘지 못하면 집합적으로 약한 도구변수라고 평가한다.

## Method of Moments을 이용한 도구변수 추정

MM은 모집단의 적률 조건이 표본에서도 그대로 성립하도록 표본 적률을 0으로 놓고 파라미터를 추정하는 직관적인 방법이다.

다음과 같은 단일 회귀 모형을 가정하자.

$$
Y_t = \beta_0 + \beta_1 X_t + u_t
$$

&#8203;$$X_t$$가 내생성을 가져 $$E[X_t u_t] \neq 0$$이지만, 타당한 도구변수 $$Z_t$$가 존재하여 다음의 모집단 직교 조건이 성립한다고 가정한다.

$$
\mathbb{E}[u_t] = 0 \quad \text{and} \quad \mathbb{E}[Z_t u_t] = 0
$$

오차항 $$u_t = Y_t - \beta_0 - \beta_1 X_t$$를 대입하면 두 개의 모집단 적률 조건을 얻는다.

$$
\begin{aligned}
\mathbb{E}[Y_t - \beta_0 - \beta_1 X_t] &= 0 \\
\mathbb{E}[Z_t (Y_t - \beta_0 - \beta_1 X_t)] &= 0
\end{aligned}
$$

이를 표본 대응치로 대체하여 0으로 둔다.

$$
\begin{aligned}
\frac{1}{n} \sum_{t=1}^n (Y_t - \hat{\beta}_0 - \hat{\beta}_1 X_t) &= 0 \\
\frac{1}{n} \sum_{t=1}^n Z_t (Y_t - \hat{\beta}_0 - \hat{\beta}_1 X_t) &= 0
\end{aligned}
$$

첫 번째 식에서 $$\hat{\beta}_0 = \bar{Y} - \hat{\beta}_1 \bar{X}$$가 도출된다. 이를 두 번째 식에 대입하여 정리한다.

$$
\begin{aligned}
\sum_{t=1}^n Z_t \left( (Y_t - \bar{Y}) - \hat{\beta}_1 (X_t - \bar{X}) \right) &= 0 \\
\sum_{t=1}^n Z_t (Y_t - \bar{Y}) &= \hat{\beta}_1 \sum_{t=1}^n Z_t (X_t - \bar{X})
\end{aligned}
$$

이를 편차 형태로 묶어주면 최종적으로 다음과 같은 IV estimator가 유도된다.

$$
\hat{\beta}_1^{IV} = \frac{\sum_{t=1}^n (Z_t - \bar{Z})(Y_t - \bar{Y})}{\sum_{t=1}^n (Z_t - \bar{Z})(X_t - \bar{X})} = \frac{\widehat{Cov}(Z, Y)}{\widehat{Cov}(Z, X)}
$$
