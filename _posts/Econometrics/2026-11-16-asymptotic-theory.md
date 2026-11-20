---
title: "[Econometrics] Asymptotic Theory"
date: 2026-11-16 17:51:00 +0900
categories: [Econometrics]
tags: [점근이론, CLT, Slutsky, OLS]
math: true
toc: true
---

## Definition

표본의 크기 $$n$$이 무한대로 커질 때, 확률변수열의 분포가 특정한 확률분포에 한없이 가까워지는 현상을 의미한다. 계량경제학에서 소표본(finite sample)의 정확한 분포를 알 수 없을 때, 대표본(large sample) 근사를 통해 추정량의 가설검정을 수행하기 위한 이론적 토대다.

**정의.** 확률변수열 $$X_1, X_2, \dots, X_n$$과 이에 대응하는 누적분포함수(CDF) $$F_n(x)$$가 있다고 하자. 어떤 확률변수 $$X$$의 누적분포함수 $$F(x)$$의 모든 연속점에서 다음이 성립할 때, 확률변수열 $$X_n$$은 $$X$$로 분포 수렴한다고 정의한다.

$$
\lim_{n \to \infty} F_n(x) = F(x)
$$

표기는 주로 $$X_n \xrightarrow{d} X$$를 사용한다.

## 유용한 정리들

분포 수렴을 증명하고 계량 모형에 적용하기 위해서는 다음의 정리들이 결합되어 사용된다.

#### 중심극한정리 (Central Limit Theorem, CLT)

유한한 분산을 가지는 i.i.d. 확률변수들의 표본평균 $$\bar{X}_n$$에 대하여, 이를 표준화한 값은 $$n$$이 커짐에 따라 표준정규분포로 분포 수렴한다.

$$
\sqrt{n}(\bar{X}_n - \mu) \xrightarrow{d} N(0, \sigma^2)
$$

#### 슬루츠키 정리 (Slutsky's Theorem)

확률 수렴과 분포 수렴을 결합하는 데 사용된다. $$X_n \xrightarrow{d} X$$이고 $$Y_n \xrightarrow{p} c$$ ($$c$$는 상수)일 때, 다음이 성립한다.

$$
\begin{aligned}
X_n + Y_n &\xrightarrow{d} X + c \\
X_n Y_n &\xrightarrow{d} cX
\end{aligned}
$$

#### 연속사상정리 (Continuous Mapping Theorem, CMT)

&#8203;$$X_n \xrightarrow{d} X$$이고 함수 $$g(\cdot)$$가 연속함수라면, $$g(X_n) \xrightarrow{d} g(X)$$가 성립한다는 정리다. 이는 확률 수렴에도 동일하게 적용된다.

## OLS 추정량의 점근적 정규분포 도출

#### OLS 추정량의 오차식 전개 및 스케일링

&#8203;$$\hat{\beta}$$의 공식에 $$Y = X\beta + u$$를 대입하여 참값 $$\beta$$와의 오차를 정리한다.

$$
\hat{\beta} = (X^\mathsf{T}X)^{-1}X^\mathsf{T}(X\beta + u) = \beta + (X^\mathsf{T}X)^{-1}X^\mathsf{T}u
$$

$$
\hat{\beta} - \beta = (X^\mathsf{T}X)^{-1}X^\mathsf{T}u
$$

극한($$n \to \infty$$)에서 분산이 0으로 수렴하지 않도록 양변에 $$\sqrt{n}$$을 곱하여 스케일링한다. 이후 확률 수렴하는 부분과 분포 수렴하는 부분으로 나누기 위해 $$n$$을 교차 배치한다.

$$
\sqrt{n}(\hat{\beta} - \beta) = \left( \frac{1}{n} X^\mathsf{T}X \right)^{-1} \left( \frac{1}{\sqrt{n}} X^\mathsf{T}u \right)
$$

#### LLN과 CMT 적용

대수의 법칙에 의해 개별 관측치 외적의 표본평균 $$\frac{1}{n} X^\mathsf{T}X$$는 모평균 행렬 $$Q$$로 확률 수렴한다.

$$
\frac{1}{n} X^\mathsf{T}X \xrightarrow{p} \mathbb{E}[x_i x_i^\mathsf{T}] \equiv Q
$$

여기에 연속사상정리를 적용하면 역행렬 연산에 대해서도 확률 수렴이 그대로 유지된다.

$$
\left( \frac{1}{n} X^\mathsf{T}X \right)^{-1} \xrightarrow{p} Q^{-1}
$$

#### CLT 적용

강외생성 가정이 성립하여 $$\mathbb{E}[x_i u_i] = 0$$이라 하자. 중심극한정리에 의해 $$\frac{1}{\sqrt{n}} \sum x_i u_i$$는 정규분포로 분포 수렴한다. 분산-공분산 행렬을 $$\Omega = \mathbb{E}[u_i^2 x_i x_i^\mathsf{T}]$$라 정의한다.

$$
\frac{1}{\sqrt{n}} X^\mathsf{T}u \xrightarrow{d} N(0, \Omega)
$$

#### Slutsky's Theorem 결합

확률 수렴하는 행렬 $$Q^{-1}$$과 분포 수렴하는 확률벡터 $$N(0, \Omega)$$의 곱에 슬루츠키 정리를 적용하여 최종 점근 분포를 도출한다.

$$
\sqrt{n}(\hat{\beta} - \beta) \xrightarrow{d} Q^{-1} \cdot N(0, \Omega) = N(0, Q^{-1} \Omega Q^{-1})
$$

가운데 $$\Omega$$를 $$Q^{-1}$$이 양옆에서 감싸는 이 형태를 샌드위치 분산(Sandwich variance)이라 부른다.

#### 오차항의 가정에 따른 점근 분산의 형태

**① 동분산성(Homoskedasticity) 가정**: 오차항의 조건부 분산이 $$\mathbb{E}[u_i^2 | x_i] = \sigma^2$$로 일정하다면 $$\Omega = \sigma^2 \mathbb{E}[x_i x_i^\mathsf{T}] = \sigma^2 Q$$가 된다. 이를 샌드위치 공식에 대입하면 약분되어 다음과 같이 단순화된다.

$$
Q^{-1} (\sigma^2 Q) Q^{-1} = \sigma^2 Q^{-1} \implies \sqrt{n}(\hat{\beta} - \beta) \xrightarrow{d} N(0, \sigma^2 Q^{-1})
$$

**② 이분산성(Heteroskedasticity) 존재**: 조건부 분산이 상수가 아니라면 $$\Omega$$를 $$\sigma^2 Q$$로 분리할 수 없으므로, $$Q^{-1} \Omega Q^{-1}$$ 형태를 그대로 유지해야 한다. 표본 데이터로 $$\Omega$$를 일치성 있게 추정하기 위해 $$\hat{\Omega} = \frac{1}{n} \sum_{i=1}^n \hat{u}_i^2 x_i x_i^\mathsf{T}$$를 대입하는 방식이 White 이분산 일치 표준오차의 도출 원리다.
