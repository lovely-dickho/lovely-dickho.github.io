---
title: "[Econometrics] Vector AutoRegression"
date: 2026-11-16 19:00:00 +0900
categories: [Econometrics]
tags: [VAR, Granger, 충격반응, 분산분해]
math: true
toc: true
---

## 교차상관분석

교차상관성은 두 시계열변수의 시차값들 사이의 상관관계를 측정하는 것이다. 두 개의 정상적시계열변수 $$X_t$$와 $$Y_t$$가 있다고 가정하자. $$k$$시차에서 $$X_t$$와 $$Y_t$$의 교차상관성은 다음과 같다.

$$
\rho_{XY}(k) = \frac{\mathbb{E}[(X_t - \mu_X)(Y_{t-k} - \mu_Y)]}{\sigma_X \sigma_Y}
\tag{1}
$$

표본 교차상관성은 다음과 같다.

$$
\hat{\rho}_{XY}(k) = \frac{\frac{1}{T-k} \sum_{t=k+1}^T (x_t - \bar{x})(y_{t-k} - \bar{y})}{\sqrt{\frac{1}{T} \sum_{t=1}^T (x_t - \bar{x})^2} \sqrt{\frac{1}{T} \sum_{t=1}^T (y_t - \bar{y})}}
\tag{2}
$$

- $$\rho_{XY}(k)$$와 $$\rho_{YX}(k)$$는 동일하지 않으므로, 첨자의 순서가 중요하다. $$\rho_{XY}(-k)$$는 $$X_t$$와 $$Y_{t+k}$$의 상관관계를 나타내므로, $$X_{t-k}$$와 $$Y_t$$의 상관관계를 의미하는 $$\rho_{YX}(k)$$와 동일하다.
- 양의 수 $$k$$에 대해 $$\rho_{XY}(k)$$가 큰 값을 갖는다는 것은 $$X_t$$와 $$Y_{t-k}$$의 상관관계가 높다는 것을 의미한다. 이는 $$Y$$의 움직임이 $$X$$를 $$k$$기간만큼 선도(lead)한다고 해석할 수 있다.

## Granger causality test

그랜저 인과성검증은 한 변수가 다른 변수의 예측에 도움이 되는지를 검증하는 방법이다.

**정의.** $$X_t$$의 시차값들이 다른 변수 $$Y_t$$를 예측하는 데 도움이 된다면, ‘$$X$$는 $$Y$$를 그랜저 인과한다’라고 말한다. 반면에 예측에 도움이 되지 않는다면 ‘$$X$$는 $$Y$$를 그랜저 인과하지 않는다’라고 말한다.

그랜저 인과성에 대한 모형은 다음과 같은 2개의 자기회귀모형이다.

$$
\begin{aligned}
Y_t &= \alpha_{10} + \alpha_{11}Y_{t-1} + \cdots + \alpha_{1p}Y_{t-p} + \beta_{11}X_{t-1} + \cdots + \beta_{1p}X_{t-p} + \varepsilon_{1t} \qquad \text{(3)} \\
Y_t &= \alpha_{20} + \alpha_{21}Y_{t-1} + \cdots + \alpha_{2p}Y_{t-p} + \varepsilon_{2t} \qquad \text{(4)}
\end{aligned}
$$

(3)의 모수 $$\beta_{11}, \cdots, \beta_{1p}$$가 0이 아니라면, $$X$$는 $$Y$$를 그랜저 인과하는 것이므로, $$X_t$$의 시차값들은 $$Y_t$$를 예측하는 데 유용하다. (4)은 모수 $$\beta_{11}, \cdots, \beta_{1p}$$가 모두 0이라는 제약하의 제약모형이다.

#### 가설 검증

그랜저 인과성검증은 $$X$$가 $$Y$$를 그랜저 인과하지 않는다는 귀무가설($$H_0 : \beta_{11} = \cdots = \beta_{1p} = 0$$)을 검증하는 것이다. (3)에서 $$X_{t-1}, \cdots, X_{t-p}$$가 독립변수로 포함됨으로써 적합도가 증가한다면, 위의 귀무가설은 기각되어야 한다.

그랜저 인과성검증을 실행하기 위해서는 OLS에 의해 제약모형인 (4)를 먼저 추정하고 $$R_R^2$$을 계산하여야 한다. 차수 $$p$$를 결정하는 데에는 아무런 제약과 규칙도 없지만 $$p$$를 충분히 큰 값으로 선택하는 것이 좋다. 만일 무제약모형인 (3)의 $$R_U^2$$이 제약모형인 (4)의 $$R_R^2$$보다 유의하게 크다면 귀무가설은 기각된다.

$$
F = \frac{(R_U^2 - R_R^2)/p}{(1 - R_U^2)/(T - 2p - 1)}
\tag{5}
$$

#### 주의

‘$$X$$가 $$Y$$를 그랜저 인과한다’라는 말은 $$X$$가 원인이고 $$Y$$가 $$X$$의 결과 또는 효과라는 것을 의미하지 않는다. 그랜저 인과성검증의 결과는 $$X$$의 선행성(precedence)과 정보내용만을 보여줄 뿐이다.

## VAR

#### 구조적(structural) 모형

피드백 효과를 갖는 두 변수 $$X_{1,t}$$와 $$X_{2,t}$$가 있다고 하면, 두 변수 사이의 관계는 다음의 2변량 시스템으로 표현된다.

$$
\begin{aligned}
X_{1,t} &= \alpha_{10} + \beta_{12} X_{2,t} + \phi_{11} X_{1,t-1} + \phi_{12} X_{2,t-1} + \varepsilon_{1,t} \qquad \text{(6)} \\
X_{2,t} &= \alpha_{20} + \beta_{21} X_{1,t} + \phi_{21} X_{1,t-1} + \phi_{22} X_{2,t-1} + \varepsilon_{2,t} \qquad \text{(7)}
\end{aligned}
$$

두 변수 $$X_{1,t}$$와 $$X_{2,t}$$는 동시적으로 결정되기 때문에 두 변수는 내생변수이다. 이 모형에서 다음이 가정된다.

**①** $$X_{1,t}$$와 $$X_{2,t}$$는 stationary하다.

**②** $$\varepsilon_{1,t}$$와 $$\varepsilon_{2,t}$$는 각각 $$\sigma_{\varepsilon_1}^2$$와 $$\sigma_{\varepsilon_2}^2$$의 분산을 갖는 백색잡음 잔차이다.

(6)와 (7)의 2변량 모형을 행렬로 표현하자.

$$
\underbrace{\begin{pmatrix} 1 & -\beta_{12} \\ -\beta_{21} & 1 \end{pmatrix}}_{\boldsymbol{\beta}} \underbrace{\begin{pmatrix} X_{1,t} \\ X_{2,t} \end{pmatrix}}_{\mathbf{X}_t} = \underbrace{\begin{pmatrix} \alpha_{10} \\ \alpha_{20} \end{pmatrix}}_{\boldsymbol{\alpha}} + \underbrace{\begin{pmatrix} \phi_{11} & \phi_{12} \\ \phi_{21} & \phi_{22} \end{pmatrix}}_{\boldsymbol{\Phi}} \underbrace{\begin{pmatrix} X_{1,t-1} \\ X_{2,t-1} \end{pmatrix}}_{\mathbf{X}_{t-1}} + \underbrace{\begin{pmatrix} \varepsilon_{1,t} \\ \varepsilon_{2,t} \end{pmatrix}}_{\boldsymbol{\varepsilon}_t}
\tag{8}
$$

&#8203;$$\mathbf{X}_t = [X_{1,t} \ X_{2,t}]^T$$는 내생변수 벡터이고, $$\boldsymbol{\varepsilon}_t = [\varepsilon_{1,t} \ \varepsilon_{2,t}]^T$$는 오차항 벡터. $$\boldsymbol{\beta}$$, $$\boldsymbol{\alpha}$$와 $$\boldsymbol{\Phi}$$는 모수 행렬이다.

2변량 모형은 $$X_{1,t}$$와 $$X_{2,t}$$ 사이의 동시적인 피드백 효과를 포함한다.

- (7)에서 $$\beta_{21}$$이 0이 아니라면, $$X_{1,t}$$에 대한 순수한 충격 $$\varepsilon_{1,t}$$는 $$X_{1,t}$$에 미치는 영향을 통해서 $$X_{2,t}$$에 대해 간접적인 동시적 효과를 미친다.
- (6)의 $$\beta_{12}$$가 0이 아니라면, $$X_{2,t}$$에 대한 순수한 충격 $$\varepsilon_{2,t}$$도 $$X_{2,t}$$에 미치는 영향을 통해서 $$X_{1,t}$$에 대해 간접적인 동시적인 효과를 미친다.

#### 축약형 모형

(8)은 $$\boldsymbol{\beta}$$의 역행렬을 양변에 곱하여 줌으로써 다음과 같은 reduced form으로 변형될 수 있다.

$$
\mathbf{X}_t = \boldsymbol{\beta}^{-1} \boldsymbol{\alpha} + \boldsymbol{\beta}^{-1} \boldsymbol{\Phi} \mathbf{X}_{t-1} + \boldsymbol{\beta}^{-1} \boldsymbol{\varepsilon}_t
$$

다음과 같이 더 단순한 형태로 나타낼 수 있다.

$$
\begin{aligned}
\underbrace{\begin{pmatrix} X_{1,t} \\ X_{2,t} \end{pmatrix}}_{\mathbf{X}_t} &=
\underbrace{\begin{pmatrix} a_{10} \\ a_{20} \end{pmatrix}}_{\mathbf{A}_0} +
\underbrace{\begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix}}_{\mathbf{A}_1}
\underbrace{\begin{pmatrix} X_{1,t-1} \\ X_{2,t-1} \end{pmatrix}}_{\mathbf{X}_{t-1}} +
\underbrace{\begin{pmatrix} e_{1,t} \\ e_{2,t} \end{pmatrix}}_{\mathbf{e}_t}
\end{aligned}
\tag{9}
$$

&#8203;$$\mathbf{A}_0 = \boldsymbol{\beta}^{-1} \boldsymbol{\alpha}$$, $$\mathbf{A}_1 = \boldsymbol{\beta}^{-1} \boldsymbol{\Phi}$$ 및 $$\mathbf{e}_t = \boldsymbol{\beta}^{-1} \boldsymbol{\varepsilon}_t$$이다.

위 식을 1차 벡터자기회귀(vector autoregression : VAR)의 표준적인 형태로서 VAR(1)모형이라고 한다.

#### VAR(p)모형

VAR(1)모형은 다음과 같은 $$k$$변량 VAR($$p$$)모형으로 일반화될 수 있다.

$$
\mathbf{X}_t = \mathbf{A}_0 + \mathbf{A}_1 \mathbf{X}_{t-1} + \cdots + \mathbf{A}_p \mathbf{X}_{t-p} + \mathbf{e}_t
\tag{10}
$$

단, $$\mathbf{X}_t = [X_{1,t} \cdots X_{k,t}]^T$$이고 $$\mathbf{e}_t = [e_{1,t} \cdots e_{k,t}]^T$$

## VAR의 추정과 식별

#### 식별 문제

일차적인 목적이 구조적 모형을 추정하는 것이라면, 표준형 VAR모형의 추정치를 가지고 구조적 모형의 모수 값을 식별해야 한다. 이때 직면하는 문제가 **식별(identification) 문제**이다.

위 모형에서 식별해야 할 모수는 8개의 계수와 2개의 잔차항의 분산으로서 총 10개이다. 축약형 모형을 추정함으로써 얻을 수 있는 추정치는 6개의 계수 추정치, 2개의 잔차항의 분산과 1개의 잔차 사이의 공분산, 합계 9개이다. 9개의 추정치를 가지고 구조적 모형의 10개의 모수를 식별해야 하므로, 과소추정(under-identification) 문제에 직면한다.

#### Cholesky decomposition

과소추정 문제를 해결하는 한 가지 대안은 구조적 모형에서 한 개의 모수값을 제약하는 것이다. 일반적인 방법은 $$\beta_{12}$$나 $$\beta_{21}$$ 중 하나를 0으로 놓음으로써, **두 변수 중 하나가 다른 하나에 동시적인 영향을 미치지 않는다고 가정**하는 것이다. 가령 (8)에서 $$X_{2,t}$$가 $$X_{1,t}$$에 동시적인 영향을 미치지 않는다고 가정하면, $$\beta_{12}$$는 0으로 제약되므로 구조적 모형은 다음과 같이 나타내진다.

$$
\begin{pmatrix} 1 & 0 \\ -\beta_{21} & 1 \end{pmatrix} \begin{pmatrix} X_{1,t} \\ X_{2,t} \end{pmatrix} = \begin{pmatrix} \alpha_{10} \\ \alpha_{20} \end{pmatrix} + \begin{pmatrix} \phi_{11} & \phi_{12} \\ \phi_{21} & \phi_{22} \end{pmatrix} \begin{pmatrix} X_{1,t-1} \\ X_{2,t-1} \end{pmatrix} + \begin{pmatrix} \varepsilon_{1,t} \\ \varepsilon_{2,t} \end{pmatrix}
\tag{11}
$$

이와 같이 모수값을 제약하는 것을 Cholesky decomposition라 한다.

## 단기 제약과 장기 제약

#### Short-run Restrictions

단기제약은 특정 구조적 충격이 특정 내생변수에 당기에 미치는 영향이 0이라고 가정하여 식별하는 방법이다. 대표적으로 Cholesky Decomposition를 이용한 Recursive 구조가 있다.

구조적 모형과 축약형 모형의 잔차항 관계 $$u_t = B_0^{-1} \varepsilon_t$$에서 행렬 $$B_0^{-1}$$를 하삼각행렬로 가정한다. 3변수 모형 $$y_t = [y_{1t}, y_{2t}, y_{3t}]^\mathsf{T}$$의 예시는 다음과 같다.

$$
\begin{bmatrix}
u_{1t} \\ u_{2t} \\ u_{3t}
\end{bmatrix}
=
\begin{bmatrix}
c_{11} & 0 & 0 \\
c_{21} & c_{22} & 0 \\
c_{31} & c_{32} & c_{33}
\end{bmatrix}
\begin{bmatrix}
\varepsilon_{1t} \\ \varepsilon_{2t} \\ \varepsilon_{3t}
\end{bmatrix}
$$

이는 $$\varepsilon_{3t}$$는 $$y_{3t}$$에만 당기에 영향을 미치며, $$\varepsilon_{1t}$$는 모든 변수에 당기에 즉각적인 영향을 미친다는 경제학적 가정(예: 통화정책의 실물경제 파급 시차)을 반영한다.

#### Long-run Restrictions

장기제약은 특정 구조적 충격이 특정 변수의 장기적인(누적된) 수준에 미치는 영향이 0이라고 가정하는 방법으로, Blanchard and Quah (1989)에 의해 널리 사용되었다.

축약형 VAR을 VMA 형태로 변환하고 구조적 충격과의 관계를 대입하면 다음과 같다.

$$
y_t = C(L) u_t = C(L) B_0^{-1} \varepsilon_t = \Theta(L) \varepsilon_t
$$

여기서 구조적 충격의 장기 누적 효과(Long-run multiplier) 행렬은 $$\Theta(1) = C(1) B_0^{-1}$$이다. 장기제약은 이 $$\Theta(1)$$ 행렬의 특정 원소를 0으로 설정하여 식별을 수행한다.

예를 들어 $$\Delta y_{1t}$$를 실질 GDP 증가율, $$y_{2t}$$를 실업률이라 하고, 충격을 공급 충격($$\varepsilon_{st}$$)과 수요 충격($$\varepsilon_{dt}$$)으로 정의할 때의 제약은 다음과 같다.

$$
\begin{bmatrix}
\Delta y_{1t} \\ \Delta y_{2t}
\end{bmatrix}
=
\begin{bmatrix}
\Theta_{11}(1) & 0 \\
\Theta_{21}(1) & \Theta_{22}(1)
\end{bmatrix}
\begin{bmatrix}
\varepsilon_{st} \\ \varepsilon_{dt}
\end{bmatrix}
$$

행렬의 우상단 원소가 0인 것은 두 번째 충격인 수요 충격($$\varepsilon_{dt}$$)이 실질 GDP($$y_{1t}$$)의 장기적인 수준에 미치는 누적 효과가 0임을 의미한다. 이는 명목 충격이 실질 변수의 장기 수준에 영향을 미치지 못한다는 장기 화폐 중립성 가설을 바탕으로 한다.

## Impulse Response Function

정상적인 AR모형이 MA($$\infty$$)모형으로 전환될 수 있는 것처럼, VAR모형도 무한한 차수의 벡터이동평균(VMA)모형으로 전환될 수 있다. (9)에서 두 변수 $$X_{1,t}$$와 $$X_{2,t}$$의 평균이 $$\mu_1$$, $$\mu_2$$일 때, 2변량 VAR(1)모형은 후방반복계산에 의해 다음의 VMA($$\infty$$)모형으로 전환될 수 있다.

$$
\begin{pmatrix} X_{1,t} \\ X_{2,t} \end{pmatrix} = \begin{pmatrix} \mu_1 \\ \mu_2 \end{pmatrix} + \sum_{l=0}^{\infty} \begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix}^l \begin{pmatrix} e_{1,t-l} \\ e_{2,t-l} \end{pmatrix}
\tag{12}
$$

(9)로부터 다음을 안다.

$$
\begin{bmatrix} e_{1,t} \\ e_{2,t} \end{bmatrix} = \begin{bmatrix} 1 & -\beta_{12} \\ -\beta_{21} & 1 \end{bmatrix}^{-1} \begin{bmatrix} \varepsilon_{1,t} \\ \varepsilon_{2,t} \end{bmatrix}
$$

이 관계를 이용하면, (12)는 구조적 모형의 잔차항인 $$\varepsilon_{1,t}$$와 $$\varepsilon_{2,t}$$의 함수로 표현될 수 있다. 그 결과 얻어지는 식을 다음과 같이 나타내자.

$$
\begin{aligned}
\begin{pmatrix} X_{1,t} \\ X_{2,t} \end{pmatrix} &= \begin{pmatrix} \mu_1 \\ \mu_2 \end{pmatrix} + \begin{pmatrix} \Psi_{11}(0) & \Psi_{12}(0) \\ \Psi_{21}(0) & \Psi_{22}(0) \end{pmatrix} \begin{pmatrix} \varepsilon_{1,t} \\ \varepsilon_{2,t} \end{pmatrix} + \begin{pmatrix} \Psi_{11}(1) & \Psi_{12}(1) \\ \Psi_{21}(1) & \Psi_{22}(1) \end{pmatrix} \begin{pmatrix} \varepsilon_{1,t-1} \\ \varepsilon_{2,t-1} \end{pmatrix} \\
&\quad + \cdots + \begin{pmatrix} \Psi_{11}(\infty) & \Psi_{12}(\infty) \\ \Psi_{21}(\infty) & \Psi_{22}(\infty) \end{pmatrix} \begin{pmatrix} \varepsilon_{1,t-\infty} \\ \varepsilon_{2,t-\infty} \end{pmatrix}
\end{aligned}
\tag{13}
$$

#### Impact Multipliers

위의 VMA($$\infty$$)모형은 충격에 대한 내생변수의 반응을 측정하는 데 실마리를 준다.

- $$X_{1,t}$$에 주어진 충격인 $$\varepsilon_{1,t}$$에 대한 $$X_{1,t}$$의 동시적인 반응은 $$\Psi_{11}(0)$$과 같다.
- $$\varepsilon_{1,t}$$에 대한 $$X_{2,t}$$의 반응은 $$\Psi_{21}(0)$$과 같고, $$\varepsilon_{2,t}$$에 대한 $$X_{1,t}$$의 반응은 $$\Psi_{12}(0)$$, $$\varepsilon_{2,t}$$에 대한 $$X_{2,t}$$의 반응은 $$\Psi_{22}(0)$$이다.
- $$\Psi_{11}(0)$$, $$\Psi_{21}(0)$$, $$\Psi_{12}(0)$$ 및 $$\Psi_{22}(0)$$는 어떤 충격에 대한 내생변수의 동시적 반응을 나타내기 때문에 이들을 충격승수라 한다.

&#8203;$$X_{j,t-l}$$에 발생했던 충격인 $$\varepsilon_{j,t-l}$$에 대한 $$X_{i,t}$$의 반응은 $$\Psi_{ij}(l)$$과 같다. 시차 $$l$$에 대한 충격승수 $$\Psi_{ij}(l)$$의 관계는 충격반응함수라 정의한다.

## Var Decomposition

(13)의 VMA($$\infty$$)모형에 포함된 정보는 (9)의 축약형 모형과 동일하기 때문에, 내생변수들의 미래 값을 예측하는 데 (13)을 이용할 수 있다.

$$
\begin{aligned}
\begin{pmatrix} X_{1,t+1} \\ X_{2,t+1} \end{pmatrix} &= \begin{pmatrix} \mu_1 \\ \mu_2 \end{pmatrix} + \begin{pmatrix} \Psi_{11}(0) & \Psi_{12}(0) \\ \Psi_{21}(0) & \Psi_{22}(0) \end{pmatrix} \begin{pmatrix} \varepsilon_{1,t+1} \\ \varepsilon_{2,t+1} \end{pmatrix} + \begin{bmatrix} \Psi_{11}(1) & \Psi_{12}(1) \\ \Psi_{21}(1) & \Psi_{22}(1) \end{bmatrix} \begin{bmatrix} \varepsilon_{1,t} \\ \varepsilon_{2,t} \end{bmatrix} + \cdots \\
&\quad + \begin{bmatrix} \Psi_{11}(\infty) & \Psi_{12}(\infty) \\ \Psi_{21}(\infty) & \Psi_{22}(\infty) \end{bmatrix} \begin{bmatrix} \varepsilon_{1,t+1-\infty} \\ \varepsilon_{2,t+1-\infty} \end{bmatrix}
\end{aligned}
$$

#### 1기 후 예측

$$
\begin{aligned}
\mathbb{E}_t \begin{pmatrix} X_{1,t+1} \\ X_{2,t+1} \end{pmatrix} &= \begin{pmatrix} \mu_1 \\ \mu_2 \end{pmatrix} + \begin{pmatrix} \Psi_{11}(1) & \Psi_{12}(1) \\ \Psi_{21}(1) & \Psi_{22}(1) \end{pmatrix} \begin{pmatrix} \varepsilon_{1,t} \\ \varepsilon_{2,t} \end{pmatrix} + \cdots + \begin{pmatrix} \Psi_{11}(\infty) & \Psi_{12}(\infty) \\ \Psi_{21}(\infty) & \Psi_{22}(\infty) \end{pmatrix} \begin{pmatrix} \varepsilon_{1,t+1-\infty} \\ \varepsilon_{2,t+1-\infty} \end{pmatrix}
\end{aligned}
$$

#### 1기후 예측오차

$$
e_t(1) = \begin{pmatrix} X_{1,t+1} \\ X_{2,t+1} \end{pmatrix} - \mathbb{E}_t \begin{pmatrix} X_{1,t+1} \\ X_{2,t+1} \end{pmatrix} = \begin{pmatrix} \Psi_{11}(0) & \Psi_{12}(0) \\ \Psi_{21}(0) & \Psi_{22}(0) \end{pmatrix} \begin{pmatrix} \varepsilon_{1,t+1} \\ \varepsilon_{2,t+1} \end{pmatrix}
$$

#### l기후 예측

동일한 방법으로, $$l$$기후 예측오차는 다음과 같이 얻어진다.

$$
\begin{aligned}
e_t(l) &= \begin{pmatrix} \Psi_{11}(0) & \Psi_{12}(0) \\ \Psi_{21}(0) & \Psi_{22}(0) \end{pmatrix} \begin{pmatrix} \varepsilon_{1,t+l} \\ \varepsilon_{2,t+l} \end{pmatrix} + \cdots + \begin{pmatrix} \Psi_{11}(l-1) & \Psi_{12}(l-1) \\ \Psi_{21}(l-1) & \Psi_{22}(l-1) \end{pmatrix} \begin{pmatrix} \varepsilon_{1,t+1} \\ \varepsilon_{2,t+1} \end{pmatrix}
\end{aligned}
\tag{14}
$$

#### l기후 예측오차분산

$$
\begin{aligned}
Var(X_{1,t+l}) &= \sigma_{\varepsilon_1}^2 \{ [\Psi_{11}(0)]^2 + \cdots + [\Psi_{11}(l-1)]^2 \} + \sigma_{\varepsilon_2}^2 \{ [\Psi_{12}(0)]^2 + \cdots + [\Psi_{12}(l-1)]^2 \}
\end{aligned}
\tag{15}
$$

#### 분산분해

(15)은 $$X_{1,t+l}$$의 예측오차분산이 두 부분으로 구성되어 있다는 것을 의미한다.

- 변수 $$X_1$$에 대한 충격인 $$\varepsilon_1$$로부터 발생한 분산
- 변수 $$X_2$$에 대한 충격인 $$\varepsilon_2$$로부터 발생한 분산

분산분해는 $$X_{1,t+l}$$의 $$l$$기후 예측오차의 전체 분산에서 이 두 가지 분산의 상대적 비중을 나타내준다.

$$
\frac{\sigma_{\varepsilon_1}^2 \{ [\Psi_{11}(0)]^2 + \cdots + [\Psi_{11}(l-1)]^2 \}}{Var(X_{1,t+l})} \quad \text{및} \quad \frac{\sigma_{\varepsilon_2}^2 \{ [\Psi_{12}(0)]^2 + \cdots + [\Psi_{12}(l-1)]^2 \}}{Var(X_{1,t+l})}
\tag{16}
$$

두 번째 부분의 비중이 모든 $$l$$에 대해 무시할 정도로 작으면, 그것은 $$\varepsilon_2$$의 충격이 모든 예측기간 동안에 걸쳐 예측오차분산을 설명하지 못한다는 것을 의미한다. 이것은 $$X_1$$이 $$\varepsilon_2$$ 및 $$X_2$$와 독립적으로 움직이고 있다는 것을 의미한다. 이와 유사한 공식이 $$X_{2,t+l}$$의 예측오차분산에 대해서도 적용될 수 있다.
