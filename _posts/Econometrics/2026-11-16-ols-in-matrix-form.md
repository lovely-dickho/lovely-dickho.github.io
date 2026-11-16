---
title: "[Econometrics] OLS in Matrix Form"
date: 2026-11-16 17:50:00 +0900
categories: [Econometrics, Notes]
tags: [OLS, GLS, FWL, 행렬, 계량경제학]
math: true
toc: true
---

## OLS Estimator

선형회귀모형은 행렬 형태로 다음과 같이 정의된다.

$$
Y = X\beta + \varepsilon
\tag{1}
$$

잔차 벡터를 $$e = Y - X\hat{\beta}$$라 할 때, 잔차제곱합(SSR)은 다음과 같이 전개된다.

$$
\begin{aligned}
SSR &= e^\mathsf{T} e = (Y - X\hat{\beta})^\mathsf{T}(Y - X\hat{\beta}) \\
&= Y^\mathsf{T} Y - \hat{\beta}^\mathsf{T} X^\mathsf{T} Y - Y^\mathsf{T} X\hat{\beta} + \hat{\beta}^\mathsf{T} X^\mathsf{T} X\hat{\beta} \\
&= Y^\mathsf{T} Y - 2\hat{\beta}^\mathsf{T} X^\mathsf{T} Y + \hat{\beta}^\mathsf{T} X^\mathsf{T} X\hat{\beta}
\end{aligned}
$$

&#8203;$$\hat{\beta}$$에 대해 일계조건을 구해 정리한다.

$$
\begin{aligned}
\frac{\partial (e^\mathsf{T} e)}{\partial \hat{\beta}} &= -2X^\mathsf{T} Y + 2X^\mathsf{T} X\hat{\beta} = 0 \\
X^\mathsf{T} X\hat{\beta} &= X^\mathsf{T} Y \quad \Longrightarrow \quad \hat{\beta} = (X^\mathsf{T} X)^{-1}X^\mathsf{T} Y
\end{aligned}
$$

#### Unbiasedness

Classical assumptions이 성립한다고 할 때, $$\hat{\beta}$$는 unbiased estimator이다.

$$
\begin{aligned}
\hat{\beta} &= (X^\mathsf{T} X)^{-1}X^\mathsf{T} Y \\
&= (X^\mathsf{T} X)^{-1}X^\mathsf{T} (X\beta + \varepsilon) \\
&= \beta + (X^\mathsf{T} X)^{-1}X^\mathsf{T} \varepsilon \\
\mathbb{E}[\hat{\beta}] &= \beta + (X^\mathsf{T} X)^{-1}X^\mathsf{T} \mathbb{E}[\varepsilon] = \beta
\end{aligned}
$$

#### Consistency

OLS 추정량의 오차항 항을 표본 크기 $$n$$으로 나누어 분해하면 다음과 같다.

$$
\begin{aligned}
\hat{\beta} &= \beta + (X^\mathsf{T} X)^{-1}X^\mathsf{T} \varepsilon \\
&= \beta + \left( \frac{1}{n} X^\mathsf{T} X \right)^{-1} \left( \frac{1}{n} X^\mathsf{T} \varepsilon \right)
\end{aligned}
$$

양변에 확률극한을 취하고 Slutsky 정리를 적용하면 다음과 같다.

$$
\begin{aligned}
\text{plim } \hat{\beta} &= \text{plim } \left[ \beta + \left( \frac{1}{n} X^\mathsf{T} X \right)^{-1} \left( \frac{1}{n} X^\mathsf{T} \varepsilon \right) \right] \\
&= \beta + \left( \text{plim } \frac{1}{n} X^\mathsf{T} X \right)^{-1} \left( \text{plim } \frac{1}{n} X^\mathsf{T} \varepsilon \right)
\end{aligned}
$$

대수의 법칙에 의해 표본 모멘트는 표본의 크기가 커짐에 따라 모집단 모멘트로 확률수렴한다.

$$
\begin{aligned}
\text{plim } \frac{1}{n} X^\mathsf{T} X &= \mathbb{E}[X_i X_i^\mathsf{T}] = \Sigma_{XX} \\
\text{plim } \frac{1}{n} X^\mathsf{T} \varepsilon &= \mathbb{E}[X_i \varepsilon_i] = 0
\end{aligned}
$$

여기서 $$\Sigma_{XX}$$는 Full Rank 가정에 의해 역행렬이 존재하는 정치행렬(Positive Definite Matrix)이며, 외생성 가정($$\mathbb{E}[\varepsilon|X]=0$$)에 의해 설명변수와 오차항의 공분산은 0이 된다.

$$
\begin{aligned}
\text{plim } \hat{\beta} &= \beta + \Sigma_{XX}^{-1} \cdot 0 = \beta
\end{aligned}
$$

따라서 OLS 추정량 $$\hat{\beta}$$는 실제 모수 $$\beta$$의 일치추정량이다.

#### Variance-covariance matrix

오차항의 분산이 homoskedasticity을 만족하고 자기상관이 없다고 가정하면($$\mathbb{E}[\varepsilon\varepsilon^\mathsf{T}] = \sigma^2 I$$), $$X$$가 non-stochastic으로 주어졌을 때 추정량 $$\hat{\beta}$$의 분산-공분산 행렬은 다음과 같이 도출된다.

$$
\begin{aligned}
Var(\hat{\beta}) &= \mathbb{E}[(\hat{\beta} - \mathbb{E}[\hat{\beta}])(\hat{\beta} - \mathbb{E}[\hat{\beta}])^\mathsf{T}] \\
&= \mathbb{E}[(X^\mathsf{T} X)^{-1}X^\mathsf{T} \varepsilon \varepsilon^\mathsf{T} X(X^\mathsf{T} X)^{-1}] \\
&= (X^\mathsf{T} X)^{-1}X^\mathsf{T} \mathbb{E}[\varepsilon \varepsilon^\mathsf{T}] X(X^\mathsf{T} X)^{-1} \\
&= (X^\mathsf{T} X)^{-1}X^\mathsf{T} (\sigma^2 I) X(X^\mathsf{T} X)^{-1} = \sigma^2 (X^\mathsf{T} X)^{-1}X^\mathsf{T} X(X^\mathsf{T} X)^{-1} \\
&= \sigma^2 (X^\mathsf{T} X)^{-1}
\end{aligned}
$$

#### Gauss-Markov Theorem

위의 고전적 가정 하에서 OLS 추정량 $$\hat{\beta}$$는 선형 불편추정량 중 가장 작은 분산을 가지는 BLUE이다.

#### 오차항 분산(σ²)의 추정과 통계적 특성

오차항의 분산 $$\sigma^2$$에 대한 최소제곱 추정량은 $$\hat{\sigma}^2 = \frac{e^\mathsf{T} e}{n-k}$$ 로 정의된다. 이때 projection matrix $$P = X(X^\mathsf{T} X)^{-1}X^\mathsf{T}$$를 이용하면 $$e = (I_n - P)Y = (I_n - P)\varepsilon$$ 이다. 대각합의 순환 성질을 이용하면 이 추정량의 불편성을 다음과 같이 증명할 수 있다.

$$
\begin{aligned}
\mathbb{E}[\hat{\sigma}^2] &= \frac{1}{n-k}\mathbb{E}[e^\mathsf{T} e] \\
&= \frac{1}{n-k}\mathbb{E}[\varepsilon^\mathsf{T}(I_n - P)\varepsilon] = \frac{1}{n-k}\mathbb{E}[\operatorname{tr}((I_n - P)\varepsilon \varepsilon^\mathsf{T})] = \frac{1}{n-k}\operatorname{tr}[(I_n - P)\mathbb{E}(\varepsilon \varepsilon^\mathsf{T})] \\
&= \frac{\sigma^2}{n-k}\operatorname{tr}(I_n - P) \\
&= \frac{\sigma^2}{n-k}(n-k) = \sigma^2
\end{aligned}
$$

여기서 $$s = \sqrt{\hat{\sigma}^2}$$ 는 회귀의 표준오차라고 부른다. 반면, $$\sigma^2$$의 MLE은 $$\hat{\sigma}^2_{MLE} = \frac{e^\mathsf{T} e}{n} = \frac{n-k}{n}\hat{\sigma}^2$$ 로 도출되며, 기댓값이 $$\mathbb{E}[\hat{\sigma}^2_{MLE}] = \frac{n-k}{n}\sigma^2 \neq \sigma^2$$ 이므로 biased 추정량이다.

#### Confidence Interval

오차항의 정규성 가정 또는 중심극한정리가 성립하는 대표본 하에서, 개별 계수 추정량 $$\hat{\beta}_j$$에 대한 $$t$$ 통계량은 다음과 같다.

$$
\begin{aligned}
se(\hat{\beta}_j) &= \sqrt{\hat{\sigma}^2 \left[ (X^\mathsf{T} X)^{-1} \right]_{jj}} \\
t &= \frac{\hat{\beta}_j - \beta_j}{se(\hat{\beta}_j)} \sim t_{n-k}
\end{aligned}
$$

유의수준 $$\alpha$$에서 모수 $$\beta_j$$에 대한 $$(1-\alpha) \times 100%$$ 신뢰구간은 다음과 같다.

$$
\begin{aligned}
\mathbb{P}\left(-t_{\alpha/2, n-k} \le \frac{\hat{\beta}_j - \beta_j}{se(\hat{\beta}_j)} \le t_{\alpha/2, n-k}\right) &= 1-\alpha \\
\implies \beta_j &\in \left[ \hat{\beta}_j - t_{\alpha/2, n-k} \cdot se(\hat{\beta}_j), \quad \hat{\beta}_j + t_{\alpha/2, n-k} \cdot se(\hat{\beta}_j) \right]
\end{aligned}
$$

대표본($$n \to \infty$$)에서는 $$t_{0.025, n-k} \approx 1.96$$이므로, 95% 신뢰구간은 $$\hat{\beta}_j \pm 1.96 \cdot se(\hat{\beta}_j)$$로 근사한다.

## Generalized Least Squares (GLS)

#### Nonspherical Errors와 OLS의 한계

선형회귀모형 (1)에서 오차항이 동분산성이나 자기상관 없음 가정을 위배하는 경우를 비구형 오차라고 한다. 이때 오차항의 기댓값과 분산-공분산 행렬은 다음과 같다.

$$
\begin{aligned}
\mathbb{E}[\varepsilon] &= 0 \\
\mathbb{E}[\varepsilon \varepsilon^\mathsf{T}] &= V \equiv \sigma^2 \Omega \neq \sigma^2 I_n
\end{aligned}
$$

비구형 오차는 heteroskedasticity ($$\Omega$$의 대각원소가 서로 다름)과 autocorrelation ($$\Omega$$의 비대각원소가 0이 아님) 등으로 발생한다.

이러한 비구형 오차 하에서도 기존의 OLS 추정량 $$\hat{\beta} = (X^\mathsf{T} X)^{-1}X^\mathsf{T} Y$$는 여전히 불편성과 일치성을 유지한다. 그러나 추정량의 분산-공분산 행렬은 다음과 같이 달라진다.

$$
\begin{aligned}
Var(\hat{\beta}) &= \mathbb{E}[(\hat{\beta} - \beta)(\hat{\beta} - \beta)^\mathsf{T}] \\
&= \mathbb{E}[(X^\mathsf{T} X)^{-1} X^\mathsf{T} \varepsilon \varepsilon^\mathsf{T} X (X^\mathsf{T} X)^{-1}] \\
&= (X^\mathsf{T} X)^{-1} X^\mathsf{T} V X (X^\mathsf{T} X)^{-1} \\
&= \sigma^2 (X^\mathsf{T} X)^{-1} X^\mathsf{T} \Omega X (X^\mathsf{T} X)^{-1} \neq \sigma^2 (X^\mathsf{T} X)^{-1}
\end{aligned}
$$

따라서 기존의 분산 공식 $$\sigma^2(X^\mathsf{T} X)^{-1}$$을 사용한 t-검정 등은 잘못된 결론을 도출할 수 있으며, OLS 추정량은 더 이상 BLUE가 아니다. 또한 오차항 분산의 통상적 추정량 $$\hat{\sigma}^2 = \frac{e^\mathsf{T} e}{n-k}$$ 역시 biased된다.

### Infeasible GLS Estimator

&#8203;$$\Omega$$ (또는 $$V$$) 행렬을 정확히 알고 있다고 가정할 때, 비구형 오차 문제를 해결하기 위한 GLS 추정량 $$\tilde{\beta}$$는 다음 두 가지 방법으로 도출 및 해석될 수 있다.

#### Weighted Least Squares

다음과 같은 $$S(b)$$를 최소화하는 해를 구한다.

$$
S(b) = (Y - Xb)^\mathsf{T} V^{-1} (Y - Xb)
$$

&#8203;$$b$$에 대해 FOC를 구하면 $$-2X^\mathsf{T} V^{-1} Y + 2X^\mathsf{T} V^{-1} X \tilde{\beta} = 0$$ 이 되며, 이를 정리하여 도출한다.

$$
\tilde{\beta} = (X^\mathsf{T} V^{-1} X)^{-1} X^\mathsf{T} V^{-1} Y = (X^\mathsf{T} \Omega^{-1} X)^{-1} X^\mathsf{T} \Omega^{-1} Y
$$

#### GLS via Transformation

&#8203;$$\Omega$$는 positive definite 대칭 행렬이므로 $$\Omega^{-1} = P P^\mathsf{T}$$를 만족하고 $$P^\mathsf{T} \Omega P = I_n$$이 성립하는 변환 행렬 $$P$$가 존재한다. 원래 모형 양변에 $$P^\mathsf{T}$$를 앞단에 곱하여 모형을 변환한다.

$$
P^\mathsf{T} Y = P^\mathsf{T} X \beta + P^\mathsf{T} \varepsilon \quad \Longrightarrow \quad Y^* = X^* \beta + \varepsilon^*
$$

이 변환된 모형에 OLS를 적용하면 위의 GLS 추정량과 완벽히 일치하는 결과를 얻는다.

$$
\begin{aligned}
\tilde{\beta} &= (X^{*\mathsf{T}} X^*)^{-1} X^{*\mathsf{T}} Y^* \\
&= (X^\mathsf{T} P P^\mathsf{T} X)^{-1} X^\mathsf{T} P P^\mathsf{T} Y \\
&= (X^\mathsf{T} \Omega^{-1} X)^{-1} X^\mathsf{T} \Omega^{-1} Y = (X^\mathsf{T} V^{-1} X)^{-1} X^\mathsf{T} V^{-1} Y
\end{aligned}
$$

이 추정량은 모수 $$\Omega$$를 사전적으로 안다고 가정하여 도출하므로 ‘infeasible GLS 추정량’이라 부른다.

### GLS 추정량의 통계적 특성

도출된 GLS 추정량 $$\tilde{\beta} = \beta + (X^\mathsf{T} V^{-1} X)^{-1} X^\mathsf{T} V^{-1} \varepsilon$$ 는 다음과 같은 통계적 특성을 지닌다.

#### Unbiasedness

$$\mathbb{E}[\tilde{\beta}] = \beta$$

#### Variance

$$
Var(\tilde{\beta}) = (X^\mathsf{T} V^{-1} X)^{-1} = \sigma^2 (X^\mathsf{T} \Omega^{-1} X)^{-1}
$$

#### Gauss-Markov Theorem (BLUE)

변환된 모형의 오차항 분산은 $$\mathbb{E}[\varepsilon^* \varepsilon^{*\mathsf{T}}] = P^\mathsf{T} (\sigma^2 \Omega) P = \sigma^2 I_n$$ 이 되어 고전적 가정을 충족한다. 따라서 GMT에 의해 GLS는 BLUE이다.

**증명.** $$\beta$$의 임의의 선형 불편추정량을 $$\bar{\beta} = A Y$$ (단, $$AX = I_k$$)라 할 때, 분산 차이는 다음과 같이 positive semi-definite를 가짐을 확인할 수 있다.

$$
\begin{aligned}
Var(\bar{\beta}) - Var(\tilde{\beta}) &= \sigma^2 \left[ A \Omega A^\mathsf{T} - (X^\mathsf{T} \Omega^{-1} X)^{-1} \right] \\
&= \sigma^2 \left[ A \Omega A^\mathsf{T} - AX (X^\mathsf{T} \Omega^{-1} X)^{-1} X^\mathsf{T} A^\mathsf{T} \right] \\
&= \sigma^2 A \Omega^{\frac{1}{2}} \left[ I_n - \Omega^{-\frac{1}{2}} X (X^\mathsf{T} \Omega^{-1} X)^{-1} X^\mathsf{T} \Omega^{-\frac{1}{2}} \right] \Omega^{\frac{1}{2}} A^\mathsf{T} \\
&= \sigma^2 A \Omega^{\frac{1}{2}} ( I_n - P_{\Omega^{-\frac{1}{2}} X} ) \Omega^{\frac{1}{2}} A^\mathsf{T} \geq 0
\end{aligned}
$$

(projection matrix $$P$$의 성질에 의해 $$I_n - P$$는 양의 준정부호 행렬이다.) ∎

#### Consistency

&#8203;$$\frac{X^\mathsf{T} V^{-1} X}{n} = \sigma^{-2} \frac{X^\mathsf{T} \Omega^{-1} X}{n}$$ 행렬이 $$n \to \infty$$ 일 때 0보다 큰 특정 유한 행렬 $$M$$으로 수렴한다는 조건이 충족되면 $$\tilde{\beta}$$는 일치추정량이 된다.

### GLS in Scalar Form

이분산성을 띠는 단순선형회귀모형 $$Y_i = \beta_0 + \beta_1 X_i + \varepsilon_i$$를 가정해 보자. 행렬 기반의 GLS 연산을 스칼라 형태로 전개하면, 이는 각 관측치에 오차항 분산의 역수($$w_i = 1/\sigma_i^2$$)를 가중치로 부여하는 WLS와 동일하다.

오차항의 분산을 상수로 만들기 위해, 모형의 양변을 관측치별 오차항의 표준편차인 $$\sigma_i$$로 나눈다.

$$
\frac{Y_i}{\sigma_i} = \beta_0 \frac{1}{\sigma_i} + \beta_1 \frac{X_i}{\sigma_i} + \frac{\varepsilon_i}{\sigma_i}
$$

이를 변환된 변수 $$Y_i^* = Y_i/\sigma_i, \; X_{0i}^* = 1/\sigma_i, \; X_{1i}^* = X_i/\sigma_i, \; \varepsilon_i^* = \varepsilon_i/\sigma_i$$ 로 치환하면 모형은 다음과 같아진다.

$$
Y_i^* = \beta_0 X_{0i}^* + \beta_1 X_{1i}^* + \varepsilon_i^*
$$

변환된 오차항 $$\varepsilon_i^*$$의 분산은 $$Var(\varepsilon_i^*) = Var\left(\frac{\varepsilon_i}{\sigma_i}\right) = \frac{1}{\sigma_i^2} Var(\varepsilon_i) = \frac{\sigma_i^2}{\sigma_i^2} = 1$$ 이 되어 동분산성을 만족하게 된다.

이 변환된 모형에 통상적인 OLS를 적용하여 잔차제곱합을 최소화하는 것은, 원래의 변수로 표현된 다음의 가중 잔차제곱합을 최소화하는 것과 같다.

$$
\begin{aligned}
SSR_{GLS} &= \sum_{i=1}^n (\varepsilon_i^*)^2 = \sum_{i=1}^n \left( \frac{Y_i - \beta_0 - \beta_1 X_i}{\sigma_i} \right)^2 \\
&= \sum_{i=1}^n w_i (Y_i - \beta_0 - \beta_1 X_i)^2
\end{aligned}
$$

이를 $$\beta_0$$와 $$\beta_1$$에 대해 각각 편미분하여 0으로 두는 FOC를 구하면 Normal Equations를 얻는다.

$$
\begin{aligned}
\sum_{i=1}^n w_i Y_i &= \hat{\beta}_0 \sum_{i=1}^n w_i + \hat{\beta}_1 \sum_{i=1}^n w_i X_i \\
\sum_{i=1}^n w_i X_i Y_i &= \hat{\beta}_0 \sum_{i=1}^n w_i X_i + \hat{\beta}_1 \sum_{i=1}^n w_i X_i^2
\end{aligned}
$$

이 연립방정식을 $$\hat{\beta}_1$$에 대해 풀면, 가중치가 반영된 최종 스칼라 추정량이 도출된다.

$$
\hat{\beta}_{1, GLS} = \frac{(\sum w_i)(\sum w_i X_i Y_i) - (\sum w_i X_i)(\sum w_i Y_i)}{(\sum w_i)(\sum w_i X_i^2) - (\sum w_i X_i)^2}
$$

이처럼 분산이 큰 관측치에는 작은 가중치를, 분산이 작은 관측치에는 큰 가중치를 부여함으로써 추정의 효율성을 극대화하는 것이 스칼라 관점에서의 GLS 핵심 원리이다.

## Projection Matrix and Residual Matrix

종속변수 $$y$$를 설명변수 행렬 $$X$$의 열공간으로 투영할 때 발생하는 잔차 벡터 $$\hat{u} = y - X\hat{\beta}$$는 $$X$$의 열공간과 직교해야 한다. 따라서 $$X$$의 모든 열벡터와 잔차 벡터의 내적은 $$0$$이 되며, 이는 다음의 직교 조건으로 표현된다.

$$
X^\mathsf{T}\hat{u} = 0
$$

위 조건식에 잔차 벡터 $$\hat{u}$$의 정의를 대입하여 전개하면 OLS 추정량 $$\hat{\beta}$$가 즉각적으로 도출된다.

$$
X^\mathsf{T}(y - X\hat{\beta}) = 0 \implies \hat{\beta} = (X^\mathsf{T}X)^{-1}X^\mathsf{T}y
$$

잔차 벡터 $$e$$는 실제 관측치와 예측치의 차이이므로, 잔차생성행렬 $$M$$으로 도출된다.

$$
\begin{aligned}
e &= Y - \hat{Y} = Y - PY = (I - P)Y = MY \\
M &= I - X(X^\mathsf{T} X)^{-1}X^\mathsf{T}
\end{aligned}
$$

잔차생성행렬 $$M$$은 대칭성($$M^\mathsf{T} = M$$)과 멱등성($$M^2 = M$$)을 가지며, 설명변수 행렬 $$X$$와 직교($$MX = 0$$)하는 성질을 지닌다.

## FWL Theorem 문제 풀이 팁

설명변수 행렬 $$X$$를 관심 있는 단일 변수 벡터 $$x_1$$과 나머지 통제변수 행렬 $$X_2$$로 분할한다. 계수 벡터 $$\beta$$ 역시 스칼라 $$\beta_1$$과 벡터 $$\beta_2$$로 나눈다.

$$
y = x_1\beta_1 + X_2\beta_2 + u
$$

나머지 변수 행렬 $$X_2$$의 열공간에 대한 직교 투영 행렬을 $$P_2 = X_2(X_2^\mathsf{T}X_2)^{-1}X_2^\mathsf{T}$$라 하고, 이에 대한 잔차 생성 행렬을 $$M_2 = I_n - P_2$$라고 정의한다. 모형의 양변에 $$M_2$$를 곱하면 다음과 같다.

$$
M_2y = M_2x_1\beta_1 + M_2X_2\beta_2 + M_2u
$$

&#8203;$$M_2$$는 $$X_2$$에 직교하므로 $$M_2X_2 = 0$$이 성립하여 식이 단순화된다.

$$
M_2y = M_2x_1\beta_1 + M_2u
$$

이 변환된 식에 최소제곱법을 적용하면 단일 계수 $$\hat{\beta}_1$$을 도출할 수 있다. $$M_2$$는 대칭 및 멱등행렬이므로 다음과 같이 정리된다.

$$
\hat{\beta}_1 = ((M_2x_1)^\mathsf{T}(M_2x_1))^{-1}(M_2x_1)^\mathsf{T}(M_2y) = (x_1^\mathsf{T}M_2x_1)^{-1}(x_1^\mathsf{T}M_2y)
$$

이는 다른 모든 통제변수들의 영향력을 $$y$$와 $$x_1$$ 양쪽에서 모두 제거한 후, $$x_1$$의 순수한 변동과 $$y$$의 순수한 변동 사이의 관계만을 측정하는 부분 효과(Partial effect)이다.

### 분산 일치 검증

#### 전체 모형에서의 분산 도출

전체 모형의 추정량 분산은 $$Var(\hat{\beta}|X) = \sigma^2(X^\mathsf{T}X)^{-1}$$이다. $$X^\mathsf{T}X$$를 블록 행렬로 나타내면 다음과 같다.

$$
X^\mathsf{T}X = \begin{bmatrix} x_1^\mathsf{T}x_1 & x_1^\mathsf{T}X_2 \\ X_2^\mathsf{T}x_1 & X_2^\mathsf{T}X_2 \end{bmatrix}
$$

&#8203;$$\hat{\beta}_1$$의 분산은 $$(X^\mathsf{T}X)^{-1}$$의 좌측 상단 요소이다. 블록 행렬의 역행렬 공식에 따라 해당 요소는 다음과 같이 도출된다.

$$
(x_1^\mathsf{T}x_1 - x_1^\mathsf{T}X_2(X_2^\mathsf{T}X_2)^{-1}X_2^\mathsf{T}x_1)^{-1} = (x_1^\mathsf{T}(I_n - P_2)x_1)^{-1} = (x_1^\mathsf{T}M_2x_1)^{-1}
$$

따라서 전체 모형에서 도출된 $$\hat{\beta}_1$$의 분산은 다음과 같다.

$$
Var(\hat{\beta}_1)_{full} = \sigma^2(x_1^\mathsf{T}M_2x_1)^{-1}
$$

#### 분할 모형에서의 분산 도출

앞서 분할 회귀에서 구한 추정량 $$\hat{\beta}_1 = (x_1^\mathsf{T}M_2x_1)^{-1}x_1^\mathsf{T}M_2y$$에 $$y = x_1\beta_1 + X_2\beta_2 + u$$를 $$M_2$$가 곱해진 형태($$M_2y = M_2x_1\beta_1 + M_2u$$)로 대입하여 실제 계수와의 차이를 구한다.

$$
\hat{\beta}_1 - \beta_1 = (x_1^\mathsf{T}M_2x_1)^{-1}x_1^\mathsf{T}M_2u
$$

이에 대한 조건부 분산을 구하면 다음과 같이 전개된다.

$$
Var(\hat{\beta}_1)_{partitioned} = \mathbb{E}\left[ \{(x_1^\mathsf{T}M_2x_1)^{-1}x_1^\mathsf{T}M_2u\} \{(x_1^\mathsf{T}M_2x_1)^{-1}x_1^\mathsf{T}M_2u\}^\mathsf{T} | X \right]
$$

$$
= (x_1^\mathsf{T}M_2x_1)^{-1} x_1^\mathsf{T}M_2 \mathbb{E}[uu^\mathsf{T}|X] M_2 x_1 (x_1^\mathsf{T}M_2x_1)^{-1}
$$

동분산성 가정 $$\mathbb{E}[uu^\mathsf{T}|X] = \sigma^2 I_n$$과 $$M_2$$의 멱등성($$M_2M_2 = M_2$$)을 적용한다.

$$
= \sigma^2 (x_1^\mathsf{T}M_2x_1)^{-1} (x_1^\mathsf{T}M_2x_1) (x_1^\mathsf{T}M_2x_1)^{-1} = \sigma^2(x_1^\mathsf{T}M_2x_1)^{-1}
$$

결과적으로 두 가지 방식으로 접근한 $$\hat{\beta}_1$$의 분산은 일치한다.

## Partitioned Regression and FWL Theorem

모형이 $$Y = X_1\beta_1 + X_2\beta_2 + \varepsilon$$으로 분할되어 있을 때, 정규방정식을 전개하면 다음과 같다.

$$
\begin{bmatrix} X_1^\mathsf{T} \\ X_2^\mathsf{T} \end{bmatrix} \begin{bmatrix} X_1 & X_2 \end{bmatrix} \begin{bmatrix} \hat{\beta}_1 \\ \hat{\beta}_2 \end{bmatrix} = \begin{bmatrix} X_1^\mathsf{T} \\ X_2^\mathsf{T} \end{bmatrix} Y
$$

행렬의 곱셈을 통해 다음 두 개의 연립방정식을 얻는다.

$$
\begin{aligned}
X_1^\mathsf{T} X_1\hat{\beta}_1 + X_1^\mathsf{T} X_2\hat{\beta}_2 &= X_1^\mathsf{T} Y \qquad \text{(2)} \\
X_2^\mathsf{T} X_1\hat{\beta}_1 + X_2^\mathsf{T} X_2\hat{\beta}_2 &= X_2^\mathsf{T} Y \qquad \text{(3)}
\end{aligned}
$$

식 (3)을 $$\hat{\beta}_2$$에 대해 풀면 다음과 같다.

$$
\hat{\beta}_2 = (X_2^\mathsf{T} X_2)^{-1}X_2^\mathsf{T} Y - (X_2^\mathsf{T} X_2)^{-1}X_2^\mathsf{T} X_1\hat{\beta}_1
\tag{4}
$$

도출된 식 (4)를 식 (2)에 대입하여 $$\hat{\beta}_2$$를 소거한다.

$$
X_1^\mathsf{T} X_1\hat{\beta}_1 + X_1^\mathsf{T} X_2 \left[ (X_2^\mathsf{T} X_2)^{-1}X_2^\mathsf{T} Y - (X_2^\mathsf{T} X_2)^{-1}X_2^\mathsf{T} X_1\hat{\beta}_1 \right] = X_1^\mathsf{T} Y
$$

이를 $$\hat{\beta}_1$$과 $$Y$$가 있는 항으로 각각 묶고, $$X_2$$에 대한 잔차생성행렬 $$M_2 = I - X_2(X_2^\mathsf{T} X_2)^{-1}X_2^\mathsf{T}$$를 적용하여 간소화한다.

$$
\begin{aligned}
X_1^\mathsf{T} \left[ I - X_2(X_2^\mathsf{T} X_2)^{-1}X_2^\mathsf{T} \right] X_1\hat{\beta}_1 &= X_1^\mathsf{T} \left[ I - X_2(X_2^\mathsf{T} X_2)^{-1}X_2^\mathsf{T} \right] Y \\
X_1^\mathsf{T} M_2 X_1 \hat{\beta}_1 &= X_1^\mathsf{T} M_2 Y \\
\therefore \hat{\beta}_1 &= (X_1^\mathsf{T} M_2 X_1)^{-1} X_1^\mathsf{T} M_2 Y
\end{aligned}
$$

이 결과는 1단계와 2단계에서 도출한 잔차 벡터 $$M_2Y$$를 잔차 행렬 $$M_2X_1$$에 대해 회귀분석하여 얻는 추정량 $$\tilde{\beta}_1$$과 수학적으로 완벽히 동일하다.

$$
\begin{aligned}
\tilde{\beta}_1 &= [(M_2X_1)^\mathsf{T}(M_2X_1)]^{-1}(M_2X_1)^\mathsf{T}(M_2Y) \\
&= [X_1^\mathsf{T} M_2^\mathsf{T} M_2 X_1]^{-1}X_1^\mathsf{T} M_2^\mathsf{T} M_2 Y \\
&= (X_1^\mathsf{T} M_2 X_1)^{-1} X_1^\mathsf{T} M_2 Y
\end{aligned}
$$

## From Matrix to Scalar Form

#### 편차 변형

설명변수가 2개인 모형($$Y_i = \beta_0 + \beta_1 X_{1i} + \beta_2 X_{2i} + \varepsilon_i$$)을 스칼라로 전개하기 위해, 대수적 연산의 편의를 도모하고자 모든 변수를 표본평균에서 차감한 편차(Deviation) 형태로 변환한다.

$$
y_i = Y_i - \bar{Y}, \quad x_{1i} = X_{1i} - \bar{X}_1, \quad x_{2i} = X_{2i} - \bar{X}_2
$$

편차 변환을 통해 상수항 $$\beta_0$$가 소거되므로, 모형은 행렬 $$y = X\beta + \varepsilon$$로 축약된다. (단, $$X$$는 $$n \times 2$$ 행렬, $$\beta$$는 $$2 \times 1$$ 벡터이다.)

$$
X = \begin{bmatrix} x_{11} & x_{21} \\ x_{12} & x_{22} \\ \vdots & \vdots \\ x_{1n} & x_{2n} \end{bmatrix}, \quad
y = \begin{bmatrix} y_1 \\ y_2 \\ \vdots \\ y_n \end{bmatrix}, \quad
\beta = \begin{bmatrix} \beta_1 \\ \beta_2 \end{bmatrix}
$$

#### 행렬의 스칼라 전개 및 역행렬 도출

투입된 편차 행렬 $$X$$와 $$y$$를 바탕으로 교차곱 행렬 $$X^\mathsf{T} X$$와 벡터 $$X^\mathsf{T} y$$를 스칼라 합계 수식으로 나타낸다.

$$
X^\mathsf{T} X = \begin{bmatrix} \sum x_{1i}^2 & \sum x_{1i}x_{2i} \\ \sum x_{1i}x_{2i} & \sum x_{2i}^2 \end{bmatrix}, \quad
X^\mathsf{T} y = \begin{bmatrix} \sum x_{1i}y_i \\ \sum x_{2i}y_i \end{bmatrix}
$$

역행렬을 구하면 다음과 같다.

$$
(X^\mathsf{T} X)^{-1} = \frac{1}{|X^\mathsf{T} X|} \begin{bmatrix} \sum x_{2i}^2 & -\sum x_{1i}x_{2i} \\ -\sum x_{1i}x_{2i} & \sum x_{1i}^2 \end{bmatrix}
$$

#### 스칼라 형태의 최종 추정량 도출

&#8203;$$\hat{\beta} = (X^\mathsf{T} X)^{-1}X^\mathsf{T} y$$ 에 앞서 구한 결과를 대입하면 다음과 같다.

$$
\begin{bmatrix} \hat{\beta}_1 \\ \hat{\beta}_2 \end{bmatrix} = \frac{1}{|X^\mathsf{T} X|} \begin{bmatrix} \sum x_{2i}^2 & -\sum x_{1i}x_{2i} \\ -\sum x_{1i}x_{2i} & \sum x_{1i}^2 \end{bmatrix} \begin{bmatrix} \sum x_{1i}y_i \\ \sum x_{2i}y_i \end{bmatrix}
$$

식을 풀어 쓰면 다음을 얻는다.

$$
\begin{aligned}
\hat{\beta}_1 &= \frac{(\sum x_{2i}^2)(\sum x_{1i}y_i) - (\sum x_{1i}x_{2i})(\sum x_{2i}y_i)}{(\sum x_{1i}^2)(\sum x_{2i}^2) - (\sum x_{1i}x_{2i})^2} \\
\hat{\beta}_2 &= \frac{(\sum x_{1i}^2)(\sum x_{2i}y_i) - (\sum x_{1i}x_{2i})(\sum x_{1i}y_i)}{(\sum x_{1i}^2)(\sum x_{2i}^2) - (\sum x_{1i}x_{2i})^2}
\end{aligned}
$$

소거되었던 상수항 $$\hat{\beta}_0$$는 OLS 회귀선이 표본평균 좌표 $$(\bar{X}_1, \bar{X}_2, \bar{Y})$$를 반드시 통과한다는 성질을 이용하여 다음과 같이 도출한다.

$$
\hat{\beta}_0 = \bar{Y} - \hat{\beta}_1 \bar{X}_1 - \hat{\beta}_2 \bar{X}_2
$$

#### 모형의 행렬 표현

다음과 같은 다중 모형을 살펴보자.

$$
y_i = \beta_0 + \beta_1 x_{1i} + \beta_2 x_{2i} + e_i \quad (i = 1, 2, \dots, n)
$$

이를 행렬 형태인 $$y = X\beta + e$$로 나타내면 각 행렬과 벡터는 다음과 같이 구성된다.

$$
\begin{bmatrix} y_1 \\ y_2 \\ \vdots \\ y_n \end{bmatrix} =
\begin{bmatrix} 1 & x_{11} & x_{21} \\ 1 & x_{12} & x_{22} \\ \vdots & \vdots & \vdots \\ 1 & x_{1n} & x_{2n} \end{bmatrix}
\begin{bmatrix} \beta_0 \\ \beta_1 \\ \beta_2 \end{bmatrix} +
\begin{bmatrix} e_1 \\ e_2 \\ \vdots \\ e_n \end{bmatrix}
$$

SSR을 최소화하기 위한 FOC를 $$\hat{\beta}$$에 대해 편미분하여 구하면 다음과 같은 Normal Equations이 도출된다.

$$
\begin{aligned}
\frac{\partial (\hat{e}^\mathsf{T} \hat{e})}{\partial \hat{\beta}} &= -2X^\mathsf{T} y + 2X^\mathsf{T} X \hat{\beta} = 0 \\
X^\mathsf{T} X \hat{\beta} &= X^\mathsf{T} y \quad \Longrightarrow \quad \hat{\beta} = (X^\mathsf{T} X)^{-1} X^\mathsf{T} y
\end{aligned}
$$

#### 2변수 모형의 구체적인 행렬 전개

앞서 정의한 행렬 $$X$$와 $$y$$를 바탕으로, 실제 연산에 대입되는 교차곱 행렬 $$X^\mathsf{T} X$$와 벡터 $$X^\mathsf{T} y$$를 스칼라 수식 전개하면 다음과 같다.

$$
\begin{aligned}
X^\mathsf{T} X &= \begin{bmatrix} 1 & 1 & \cdots & 1 \\ x_{11} & x_{12} & \cdots & x_{1n} \\ x_{21} & x_{22} & \cdots & x_{2n} \end{bmatrix}
\begin{bmatrix} 1 & x_{11} & x_{21} \\ 1 & x_{12} & x_{22} \\ \vdots & \vdots & \vdots \\ 1 & x_{1n} & x_{2n} \end{bmatrix} = \begin{bmatrix} n & \sum x_{1i} & \sum x_{2i} \\ \sum x_{1i} & \sum x_{1i}^2 & \sum x_{1i}x_{2i} \\ \sum x_{2i} & \sum x_{1i}x_{2i} & \sum x_{2i}^2 \end{bmatrix} \\
X^\mathsf{T} y &= \begin{bmatrix} 1 & 1 & \cdots & 1 \\ x_{11} & x_{12} & \cdots & x_{1n} \\ x_{21} & x_{22} & \cdots & x_{2n} \end{bmatrix}
\begin{bmatrix} y_1 \\ y_2 \\ \vdots \\ y_n \end{bmatrix} = \begin{bmatrix} \sum y_i \\ \sum x_{1i}y_i \\ \sum x_{2i}y_i \end{bmatrix}
\end{aligned}
$$

## Effect of Variable Scaling

변수의 측정 단위(Scale)를 변경할 때 추정량에 미치는 영향을 행렬 연산을 통해 확인할 수 있다. 종속변수 $$Y$$에 스칼라 상수 $$c$$를 곱하고, 각 설명변수 $$X$$에 특정 상수들이 대각선에 위치한 대각행렬 $$D$$를 곱하여 변환한다고 가정하자 ($$Y^* = cY, \; X^* = XD$$).

변환된 데이터에 기반한 새로운 추정량 $$\hat{\beta}^*$$은 다음과 같이 전개된다.

$$
\begin{aligned}
\hat{\beta}^* &= (X^{*\mathsf{T}} X^*)^{-1} X^{*\mathsf{T}} Y^* \\
&= ((XD)^\mathsf{T} (XD))^{-1} (XD)^\mathsf{T} (cY) \\
&= (D^\mathsf{T} X^\mathsf{T} X D)^{-1} D^\mathsf{T} X^\mathsf{T} (cY) \\
&= D^{-1} (X^\mathsf{T} X)^{-1} (D^\mathsf{T})^{-1} D^\mathsf{T} X^\mathsf{T} cY \\
&= c D^{-1} (X^\mathsf{T} X)^{-1} X^\mathsf{T} Y \\
&= c D^{-1} \hat{\beta}
\end{aligned}
$$

즉, 종속변수 $$Y$$에 $$c$$배 스케일링을 적용하면 추정된 회귀계수 벡터 전체도 $$c$$배가 되며, 설명변수 행렬 $$X$$의 $$k$$번째 변수에 상수 $$d_k$$를 곱하면 해당 변수에 대한 추정 계수 $$\hat{\beta}_k$$는 $$1/d_k$$배로 반비례하여 변화한다.

#### Scalar Form

변수 스케일링의 효과를 직관적으로 이해하기 위해, 설명변수가 1개인 단순선형회귀모형($$Y_i = \beta_0 + \beta_1 X_i + \varepsilon_i$$)을 스칼라 수식으로 전개해 볼 수 있다.

종속변수에 스칼라 상수 $$c$$를 곱하고($$Y_i^* = cY_i$$), 설명변수에 상수 $$d$$를 곱하여 측정 단위를 변환($$X_i^* = dX_i$$)한다고 가정하자. 변환된 변수들의 표본평균 역시 $$\bar{Y}^* = c\bar{Y}$$, $$\bar{X}^* = d\bar{X}$$가 된다.

편차의 교차곱 연산을 통해 새로운 기울기 추정량 $$\hat{\beta}_1^*$$을 전개하면 다음과 같다.

$$
\begin{aligned}
\hat{\beta}_1^* &= \frac{\sum (X_i^* - \bar{X}^*)(Y_i^* - \bar{Y}^*)}{\sum (X_i^* - \bar{X}^*)^2} \\
&= \frac{\sum (dX_i - d\bar{X})(cY_i - c\bar{Y})}{\sum (dX_i - d\bar{X})^2} \\
&= \frac{cd \sum (X_i - \bar{X})(Y_i - \bar{Y})}{d^2 \sum (X_i - \bar{X})^2} \\
&= \frac{c}{d} \left[ \frac{\sum (X_i - \bar{X})(Y_i - \bar{Y})}{\sum (X_i - \bar{X})^2} \right] \\
&= \frac{c}{d} \hat{\beta}_1
\end{aligned}
$$

새로운 절편 추정량 $$\hat{\beta}_0^*$$ 역시 OLS 회귀선이 표본평균 좌표를 반드시 통과한다는 성질을 이용하여 전개한다.

$$
\begin{aligned}
\hat{\beta}_0^* &= \bar{Y}^* - \hat{\beta}_1^* \bar{X}^* \\
&= c\bar{Y} - \left( \frac{c}{d} \hat{\beta}_1 \right) (d\bar{X}) \\
&= c\bar{Y} - c\hat{\beta}_1 \bar{X} \\
&= c(\bar{Y} - \hat{\beta}_1 \bar{X}) \\
&= c\hat{\beta}_0
\end{aligned}
$$

행렬 연산의 결론과 동일하게, 종속변수에 $$c$$를 곱하면 절편을 포함한 모든 추정량($$\hat{\beta}_0, \hat{\beta}_1$$)이 $$c$$배가 되며, 특정 설명변수에 $$d$$를 곱하면 해당 설명변수의 기울기 추정량($$\hat{\beta}_1$$)만 $$1/d$$배로 반비례하여 변화함을 스칼라 전개를 통해서도 명확히 확인할 수 있다.

## Variance Decomposition

OLS의 투영 개념을 바탕으로, 종속변수 $$y$$의 총 변동을 모형이 설명하는 변동과 설명하지 못하는 변동으로 직교 분해(orthogonal decomposition)할 수 있다.

모형 $$y = \hat{y} + \hat{u}$$의 양변에서 표본 평균 $$\bar{y}$$를 빼서 편차 형태로 만들면 다음과 같다.

$$
(y - \bar{y}i) = (\hat{y} - \bar{y}i) + \hat{u}
$$

이 식의 양변을 제곱하여 총 변동을 전개하면 교차항이 발생한다.

$$
(y - \bar{y}i)^\mathsf{T}(y - \bar{y}i) = (\hat{y} - \bar{y}i)^\mathsf{T}(\hat{y} - \bar{y}i) + \hat{u}^\mathsf{T}\hat{u} + 2(\hat{y} - \bar{y}i)^\mathsf{T}\hat{u}
$$

projection matrix의 기하학적 성질에 의해 잔차 벡터 $$\hat{u}$$는 예측치 벡터 $$\hat{y}$$가 속한 열공간 $$C(X)$$와 직교하므로, 교차항 $$2(\hat{y} - \bar{y}i)^\mathsf{T}\hat{u}$$는 $$0$$이 된다. 따라서 총 변동은 다음과 같이 완전히 두 개의 독립적인 분산 성분으로 분해된다.

$$
SST = ESS + SSR
$$

이러한 분산 분해는 결정계수 $$R^2$$ 산출 및 모형의 전반적인 유의성을 검정하는 F-통계량을 도출하는 수리적 기반이 된다.

## R² 도출

&#8203;$$y$$를 $$n \times 1$$ 종속변수 벡터, $$X$$를 $$n \times k$$ 설명변수 행렬, $$\hat{\beta} = (X^\mathsf{T}X)^{-1}X^\mathsf{T}y$$를 최소제곱추정량이라고 할 때, 표본 평균을 중심으로 한 편차 변동을 구하기 위해 편차 행렬(centering matrix) $$M^0$$를 다음과 같이 정의한다.

$$
M^0 = I_n - \frac{1}{n}ii^\mathsf{T}
$$

여기서 $$i$$는 모든 원소가 1인 $$n \times 1$$ 벡터이다. $$M^0$$는 대칭 및 멱등행렬의 특성을 가진다.

총제곱합($$SST$$)은 편차 행렬을 사용하여 다음과 같이 나타낸다.

$$
SST = \sum (y_i - \bar{y})^2 = y^\mathsf{T}M^0y
$$

추정된 잔차 벡터 $$\hat{u} = y - X\hat{\beta}$$에 대해 투영 행렬 $$P = X(X^\mathsf{T}X)^{-1}X^\mathsf{T}$$와 잔차 생성 행렬 $$M = I_n - P$$를 도입하면 $$\hat{u} = My$$가 성립한다. 잔차제곱합($$SSR$$)은 다음과 같이 정리된다.

$$
SSR = \hat{u}^\mathsf{T}\hat{u} = (My)^\mathsf{T}(My) = y^\mathsf{T}M^\mathsf{T}My = y^\mathsf{T}My
$$

결정계수 $$R^2$$는 총 변동 중 모형에 의해 설명되는 변동의 비율이므로 행렬 연산을 통해 다음과 같이 도출된다.

$$
R^2 = 1 - \frac{SSR}{SST} = 1 - \frac{y^\mathsf{T}My}{y^\mathsf{T}M^0y}
$$
