---
title: "[Econometrics] Stationary Time Series"
date: 2026-11-16 18:59:00 +0900
categories: [Econometrics, Notes]
tags: [정상성, AR, MA, ARMA, Yule-Walker, Wold]
math: true
toc: true
---

## Stationarity

정상성은 확률과정(Stochastic Process)의 확률적 특성이 시간에 따라 변하지 않는 성질을 의미한다.
정상성은 크게 강정상성(Strict Stationarity)과 약정상성(Weak Stationarity)으로 구분한다.

강정상성은 현실적으로 검증하기 매우 어렵기 때문에, 경제학에서는 주로 1차 및 2차 Moments에만 제한을 두는 약정상성을 사용한다.
다음 세 가지 조건을 모두 만족하는 확률과정 $$\{Y_t\}$$ 를 약정상과정이라고 한다.

#### 평균 (E[Yₜ] = μ)

모든 $$t$$에 대하여 평균이 일정하다.

#### 분산 (Var(Yₜ) = γ₀ < ∞)

모든 $$t$$에 대하여 분산이 유한하고 일정하다.

#### 자기공분산 (Cov(Yₜ, Yₜ₋ₖ) = γₖ)

자기공분산이 시점 $$t$$가 아닌 오직 시차 $$k$$에만 의존한다.

가우시안 확률과정(Gaussian Process)의 경우, 약정상성을 만족하면 강정상성도 만족하게 된다.
강정상성은 간단히, 확률과정의 모든 Moment뿐만 아니라 분포 전체가 시간에 대해 불변임을 의미한다.

## Wold Representation Theorem

어떤 시계열 데이터가 주어졌을 때, WRT를 통해 해당 모형을 ‘선형 모형’으로 표현할 수 있다.

**정리 (Wold Representation).** 기댓값이 0인 임의의 약정상과정 $$\{Y_t\}$$ 는 다음 식 (1)과 같이 두 가지 직교하는(orthogonal) 부분으로 유일하게 분해될 수 있다.

$$
Y_t = \underbrace{\sum_{j=0}^{\infty} \psi_j \varepsilon_{t-j}}_{\text{Stochastic component}} + \underbrace{V_t}_{\text{Deterministic component}}
\tag{1}
$$

#### Stochastic component

&#8203;$$\varepsilon_t$$는 White Noise 과정으로, $$\mathbb{E}[\varepsilon_t] = 0$$, $$\mathbb{E}[\varepsilon_t^2] = \sigma^2$$, $$\mathbb{E}[\varepsilon_t \varepsilon_s] = 0 \ (t \neq s)$$ 이다.
식 (1)의 계수 $$\psi_j$$ 는 $$\psi_0 = 1$$ 이며, $$\sum_{j=0}^{\infty} \psi_j^2 < \infty$$을 만족한다.

#### Deterministic component

&#8203;$$V_t$$ 는 과거의 값들만으로 오차 없이 완벽하게 예측 가능한 과정이다.
모든 $$t, s$$에 대해 $$\mathbb{E}[\varepsilon_t V_s] = 0$$ 이다.

거시경제학에서 $$V_t$$는 보통 상수항(평균)이나 결정론적 추세로 간주하여 Detrending할 수 있다.
따라서 우리가 분석하는 모든 정상 시계열은 무한 이동평균 과정, 즉 MA($$\infty$$) 형태로 표현할 수 있다.
이는 복잡한 AR($$p$$)나 ARMA($$p, q$$) 모형이라도 정상성만 만족한다면, 과거 충격($$\varepsilon$$)들의 선형 결합으로 전개하여 분석할 수 있음을 의미한다.

## AR (Autoregressive) 모형

AR($$p$$) 모형은 특성 방정식의 모든 근이 단위원 밖에 존재하여 정상성을 만족한다고 가정한다.

$$
Y_t = c + \phi_1 Y_{t-1} + \cdots + \phi_p Y_{t-p} + \varepsilon_t
\tag{2}
$$

식 (2)의 무조건부 기댓값 $$\mu$$는 $$\mu = \frac{c}{1 - \phi_1 - \cdots - \phi_p}$$ 이다.
평균이 차감된 변수 $$\tilde{Y}_t = Y_t - \mu$$를 도입하여 Yule-Walker 방정식을 구성한다.

### Yule-Walker Equation

차수가 $$p$$인 자기회귀모형 AR($$p$$)를 평균이 0으로 조정된 시계열 $$Y_t$$에 대해 정의하면 다음과 같다.

$$
Y_t = \phi_1 Y_{t-1} + \phi_2 Y_{t-2} + \cdots + \phi_p Y_{t-p} + \varepsilon_t = \sum_{i=1}^p \phi_i Y_{t-i} + \varepsilon_t
\tag{3}
$$

#### k = 0 일 때: 분산 (γ₀)

&#8203;$$Y_t$$와 자기 자신인 $$Y_t$$의 공분산(분산)을 구하기 위해 식 (3)의 양변에 $$Cov(\cdot, Y_t)$$를 취한다.

$$
\begin{aligned}
Cov(Y_t, Y_t) &= Cov\left( \sum_{i=1}^p \phi_i Y_{t-i} + \varepsilon_t, Y_t \right) \\
&= \sum_{i=1}^p \phi_i Cov(Y_{t-i}, Y_t) + Cov(\varepsilon_t, Y_t)
\end{aligned}
$$

여기서 $$Cov(Y_t, Y_t) = \gamma_0$$, 시차 $$i$$에 대해 $$Cov(Y_{t-i}, Y_t) = \gamma_i$$이다.
현재 오차항은 현재 $$Y_t$$에만 영향을 미치므로 $$Cov(\varepsilon_t, Y_t) = \sigma^2$$가 되어 다음 식 (4)이 도출된다.

$$
\gamma_0 = \phi_1 \gamma_1 + \phi_2 \gamma_2 + \cdots + \phi_p \gamma_p + \sigma^2
\tag{4}
$$

#### k > 0 일 때: 자기공분산 (γₖ)

시차가 $$k$$만큼 나는 $$Y_{t-k}$$와의 공분산을 구하기 위해 식 (3)의 양변에 $$Cov(\cdot, Y_{t-k})$$를 취한다.

$$
\begin{aligned}
Cov(Y_t, Y_{t-k}) &= Cov\left( \sum_{i=1}^p \phi_i Y_{t-i} + \varepsilon_t, Y_{t-k} \right) \\
&= \sum_{i=1}^p \phi_i Cov(Y_{t-i}, Y_{t-k}) + Cov(\varepsilon_t, Y_{t-k})
\end{aligned}
$$

&#8203;$$Cov(Y_t, Y_{t-k}) = \gamma_k$$이며, 공분산의 대칭성에 의해 $$Cov(Y_{t-i}, Y_{t-k}) = \gamma_{k-i}$$이다.
과거 관측치와 현재 오차항은 독립이므로 $$Cov(\varepsilon_t, Y_{t-k}) = 0$$이 되어 다음의 차분방정식 (5)이 성립한다.

$$
\gamma_k = \phi_1 \gamma_{k-1} + \phi_2 \gamma_{k-2} + \cdots + \phi_p \gamma_{k-p} \quad (k \ge 1)
\tag{5}
$$

#### 자기상관계수 (ACF)와 Yule-Walker 방정식

도출된 식 (5)의 양변을 $$\gamma_0$$로 나누면 자기상관계수 $$\rho_k = \frac{\gamma_k}{\gamma_0}$$에 대한 식 (6)가 얻어진다.

$$
\rho_k = \phi_1 \rho_{k-1} + \phi_2 \rho_{k-2} + \cdots + \phi_p \rho_{k-p} \quad (k \ge 1)
\tag{6}
$$

이 식에 $$k = 1, 2, \dots, p$$를 순차적으로 대입하여 전개한 $$p$$개의 선형 연립방정식을 행렬 형태로 나타낸 것이 Yule-Walker 방정식 (7)이다. ($$\rho_0 = 1, \rho_{-j} = \rho_j$$)

$$
\begin{bmatrix} \rho_1 \\ \rho_2 \\ \vdots \\ \rho_p \end{bmatrix} = \begin{bmatrix} 1 & \rho_1 & \cdots & \rho_{p-1} \\ \rho_1 & 1 & \cdots & \rho_{p-2} \\ \vdots & \vdots & \ddots & \vdots \\ \rho_{p-1} & \rho_{p-2} & \cdots & 1 \end{bmatrix} \begin{bmatrix} \phi_1 \\ \phi_2 \\ \vdots \\ \phi_p \end{bmatrix}
\tag{7}
$$

추가로, $$k=0$$일 때 구했던 식 (4)에 $$\gamma_i = \gamma_0 \rho_i$$를 대입하여 정리하면 AR($$p$$) 모형의 분산 (8)을 도출할 수 있다.

$$
\begin{aligned}
\gamma_0 &= \gamma_0 (\phi_1 \rho_1 + \phi_2 \rho_2 + \cdots + \phi_p \rho_p) + \sigma^2  \\
Var(Y_t) &= \gamma_0 = \frac{\sigma^2}{1 - \sum_{i=1}^p \phi_i \rho_i} \qquad \text{(8)}
\end{aligned}
$$

### AR 모형에서의 교차곱 처리

Yule-Walker 방정식을 세우기 위해 양변에 $$\tilde{Y}_{t-k}$$를 곱하고 기댓값을 취할 때 발생하는 교차곱 $$\mathbb{E}[\varepsilon_t \tilde{Y}_{t-k}]$$는 $$k$$에 따라 다르게 처리된다.

#### k ≥ 1 (과거 변수와의 교차곱)

&#8203;$$\tilde{Y}_{t-k}$$는 정리 (1)의 WRT에 의해 과거 오차항들($$\varepsilon_{t-k}, \varepsilon_{t-k-1}, \dots$$)의 선형결합이다.
현재의 충격 $$\varepsilon_t$$는 과거 오차항들과 직교하므로 교차곱은 완전히 소거된다.

$$
\mathbb{E}[\varepsilon_t \tilde{Y}_{t-k}] = 0
$$

#### k = 0 (현재 변수와의 교차곱)

분산을 구하기 위해 $$\tilde{Y}_t$$를 곱할 때, $$\tilde{Y}_t$$를 AR 원형으로 치환하여 전개한다.

$$
\begin{aligned}
\mathbb{E}[\varepsilon_t \tilde{Y}_t] &= \mathbb{E}[\varepsilon_t (\phi_1 \tilde{Y}_{t-1} + \cdots + \varepsilon_t)] \\
&= \phi_1 \mathbb{E}[\varepsilon_t \tilde{Y}_{t-1}] + \cdots + \mathbb{E}[\varepsilon_t^2] = \mathbb{E}[\varepsilon_t^2] = \sigma^2
\end{aligned}
$$

## MA(q) 모형

차수가 $$q$$인 이동평균모형 MA($$q$$)는 현재의 시계열 값이 평균과 현재 및 과거 $$q$$기까지의 백색잡음 오차항들의 선형결합으로 표현되는 모형이다.
수식 전개의 편의를 위해 $$\theta_0 = 1$$로 정의한다.

$$
Y_t = \mu + \varepsilon_t + \theta_1 \varepsilon_{t-1} + \cdots + \theta_q \varepsilon_{t-q} = \mu + \sum_{i=0}^q \theta_i \varepsilon_{t-i}
\tag{9}
$$

#### 무조건부 기댓값

오차항 $$\varepsilon_t$$는 백색잡음이므로 모든 시점에서 기댓값이 0이다. $$\mathbb{E}$$의 선형성에 의해 $$\mathbb{E}[Y_t] = \mathbb{E}\left[\mu + \sum_{i=0}^q \theta_i \varepsilon_{t-i}\right] = \mu + \sum_{i=0}^q \theta_i \mathbb{E}[\varepsilon_{t-i}] = \mu$$가 성립한다.

#### 분산 및 자기공분산 (γₖ)

모든 교차곱 처리의 근본적인 근거는 백색잡음 가정이다.
서로 다른 시점의 오차항 간 교차곱의 기댓값은 모두 0이 되고($$\mathbb{E}[\varepsilon_t \varepsilon_s] = 0, t \neq s$$), 시점이 일치하는 오차항의 제곱에 대한 기댓값만 분산 $$\sigma^2$$으로 살아남는다.

##### 분산 (k=0)

제곱을 전개할 때 발생하는 모든 교차곱이 소거되어 식 (10)이 도출된다.

$$
\gamma_0 = Var(Y_t) = \mathbb{E}\left[\left(\sum_{i=0}^q \theta_i \varepsilon_{t-i}\right)^2\right] = \sigma^2 \sum_{i=0}^q \theta_i^2
\tag{10}
$$

##### 자기공분산 (1 ≤ k ≤ q)

두 다항식을 곱할 때 시점이 일치하는 오차항($$i = j+k$$)끼리 짝지어질 때만 $$\sigma^2$$ 값을 갖는다.

$$
\begin{aligned}
\gamma_k &= \mathbb{E}[(Y_t - \mu)(Y_{t-k} - \mu)]  \\
&= \mathbb{E}\left[\left(\sum_{i=0}^q \theta_i \varepsilon_{t-i}\right)\left(\sum_{j=0}^q \theta_j \varepsilon_{t-k-j}\right)\right] = \sigma^2 \sum_{i=0}^{q-k} \theta_{i+k}\theta_i \qquad \text{(11)}
\end{aligned}
$$

##### 절단성 (k > q)

시차 $$k$$가 $$q$$를 초과하면, 겹치는 시점의 오차항이 단 하나도 존재하지 않게 되므로 교차곱의 기댓값이 완벽하게 0이 된다.

$$
\gamma_k = 0 \quad (k > q)
$$

#### 자기상관계수 (ACF, ρₖ)

자기상관계수는 $$\rho_k = \frac{\gamma_k}{\gamma_0}$$로 정의되므로 식 (10)와 (11)에 의해 다음과 같이 일반화된다.

$$
\begin{aligned}
\rho_k &= \frac{\sum_{i=0}^{q-k} \theta_{i+k}\theta_i}{\sum_{i=0}^q \theta_i^2} \quad (1 \le k \le q) \qquad \text{(12)} \\
\rho_k &= 0 \quad (k > q)
\end{aligned}
$$

MA($$q$$) 모형은 이처럼 시차 $$q$$ 이후로 자기상관계수가 0으로 완전히 떨어지는 **절단성(Cut-off)**을 가지며, 이는 실제 데이터 분석 시 모형의 차수 $$q$$를 식별하는 핵심 기준이 된다.

### MLE

시계열 데이터는 상호 독립이 아니므로 확률의 chain rule을 이용한 prediction error decomposition 형태로 결합 확률밀도함수를 전개한다.

$$
f(y_1, y_2, \cdots, y_T; \boldsymbol{\theta}) = f(y_1, \cdots, y_p; \boldsymbol{\theta}) \prod_{t=p+1}^T f(y_t | y_{t-1}, \cdots, y_1; \boldsymbol{\theta})
\tag{13}
$$

#### AR 모형의 조건부 MLE

초기 관측치 $$(y_1, \dots, y_p)$$ 가 상수로 주어졌다고 취급하면 식 (13)의 결합 밀도함수 앞부분을 무시할 수 있다.
가우시안 오차를 가정할 때, AR(3)의 조건부 로그 우도함수 $$\mathcal{L}_c$$ 는 다음과 같다.

$$
\mathcal{L}_c(\boldsymbol{\phi}, \sigma^2) = -\frac{T-3}{2}\ln(2\pi) - \frac{T-3}{2}\ln(\sigma^2) - \frac{1}{2\sigma^2} \sum_{t=4}^T (y_t - \phi_1 y_{t-1} - \phi_2 y_{t-2} - \phi_3 y_{t-3})^2
\tag{14}
$$

이를 극대화하는 것은 RSS를 최소화하는 것과 동일하므로, AR 모형에서 조건부 MLE는 OLS 추정량과 일치한다.

#### AR 모형의 정확한 MLE

표본 크기가 작을 때 초기 관측치의 unconditional distribution를 우도함수에 포함시키는 방식이다.
AR(1) $$y_t = \phi y_{t-1} + \varepsilon_t$$ 의 경우, $$y_1 \sim \mathcal{N}\left(0, \frac{\sigma^2}{1-\phi^2}\right)$$ 임을 반영한 정확한 로그 우도함수 $$\mathcal{L}_E$$ 는 식 (15)과 같다.

$$
\mathcal{L}_E(\phi, \sigma^2) = \underbrace{ \frac{1}{2}\ln(1-\phi^2) - \frac{1-\phi^2}{2\sigma^2}y_1^2 }_{f(y_1) \text{에서 도출된 추가항}} - \frac{T-1}{2}\ln(2\pi) - \frac{T}{2}\ln(\sigma^2) - \frac{1}{2\sigma^2} \sum_{t=2}^T (y_t - \phi y_{t-1})^2
\tag{15}
$$

앞에 추가된 비선형 항 $$\frac{1}{2}\ln(1-\phi^2)$$ 때문에 OLS로 풀 수 없으며 수치적 최적화가 필요하다.
(단, $$T \to \infty$$ 이면 OLS와 점근적 동치이다.)

#### ARMA 모형의 MLE와 잔차의 재귀적 연산

ARMA(1,1) 모형 $$y_t = \phi y_{t-1} + \varepsilon_t + \theta \varepsilon_{t-1}$$ 과 같이 MA 항이 포함되면 오차항 $$\varepsilon_{t-1}$$ 이 관측되지 않아 $$\mathbf{X}$$ 를 구성할 수 없다.
조건부 MLE를 위해 $$\varepsilon_0 = 0$$ 으로 가정하고, 잔차 $$\hat{\varepsilon}_t$$ 를 다음과 같이 재귀적으로 계산한다.

$$
\begin{aligned}
\hat{\varepsilon}_1 &= y_1 \\
\hat{\varepsilon}_2 &= y_2 - \phi y_1 - \theta \hat{\varepsilon}_1 \\
\hat{\varepsilon}_3 &= y_3 - \phi y_2 - \theta \hat{\varepsilon}_2 \\
&\vdots
\end{aligned}
$$

이 재귀적 잔차들을 모아 조건부 우도함수의 $$\sum \hat{\varepsilon}_t^2$$ 에 대입하여 최적화하며, 이는 본질적으로 비선형 최소제곱법(NLS)의 구조를 갖는다.

## 예측

분석의 편의를 위해 다음과 같은 ARMA(1,1) 모형 (16)을 가정한다. 백색잡음 가정에 의해 $$Var(\varepsilon_t) = \sigma^2$$이며, 서로 다른 시점의 오차항 간 공분산은 0이다.

$$
Y_t = c + \phi Y_{t-1} + \varepsilon_t + \theta \varepsilon_{t-1}
\tag{16}
$$

### 조건부 기댓값 연산의 기본 원리

시계열에서 미래를 forecasting하는 것은 현재까지 주어진 정보를 바탕으로 미래 변수의 조건부 기댓값을 구하는 과정이다.
현재 시점 $$t$$까지 관측된 모든 정보의 집합을 $$\Omega_t = \{Y_t, Y_{t-1}, \cdots, \varepsilon_t, \varepsilon_{t-1}, \cdots\}$$라 하자.

#### 과거 및 현재의 변수 (상수 취급)

시점 $$t$$ 이하의 모든 $$Y_{t-j}$$와 오차항 $$\varepsilon_{t-j}$$ ($$j \ge 0$$)는 이미 실현되어 관측된 값이므로 불확실성이 없다.
따라서 조건부 기댓값을 취하면 자기 자신이 그대로 나온다.

$$
\begin{aligned}
\mathbb{E}[Y_{t-j} | \Omega_t] &= Y_{t-j} \\
\mathbb{E}[\varepsilon_{t-j} | \Omega_t] &= \varepsilon_{t-j}
\end{aligned}
$$

#### 미래의 오차항 (0으로 소거)

미래에 발생할 충격 $$\varepsilon_{t+h}$$ ($$h \ge 1$$)는 현재 시점의 정보 집합 $$\Omega_t$$와 완전히 독립적인 백색잡음이다.
따라서 현재 시점에서는 미래의 충격을 전혀 예측할 수 없으므로 기댓값은 0이 된다.

$$
\mathbb{E}[\varepsilon_{t+h} | \Omega_t] = 0
$$

#### 미래의 변수 (예측치로 치환)

미래의 시계열 값 $$Y_{t+h}$$에 대한 조건부 기댓값 자체가 우리가 구하고자 하는 $$h$$기 앞 예측치(Forecast) $$\hat{Y}_{t+h|t}$$이다.

$$
\mathbb{E}[Y_{t+h} | \Omega_t] = \hat{Y}_{t+h|t}
$$

### 1기 앞 예측

1기 앞의 실제 데이터 생성 식은 (17)과 같다.

$$
Y_{t+1} = c + \phi Y_t + \varepsilon_{t+1} + \theta \varepsilon_t
\tag{17}
$$

#### 예측치 (Ŷ(t+1|t))

식 (17)의 양변에 조건부 기댓값 $$\mathbb{E}[\cdot | \Omega_t]$$를 취하면 미래의 충격 $$\varepsilon_{t+1}$$은 0으로 소거되어 식 (18)이 도출된다.

$$
\hat{Y}_{t+1|t} = c + \phi Y_t + \theta \varepsilon_t
\tag{18}
$$

#### 예측오차 (e(t+1|t))

실제값에서 예측치를 뺀 결과는 식 (19)과 같다.

$$
e_{t+1|t} = Y_{t+1} - \hat{Y}_{t+1|t} = \varepsilon_{t+1}
\tag{19}
$$

#### 예측오차의 분산

예측오차가 미래의 오차항 하나로만 구성되므로 다음과 같다.

$$
Var(e_{t+1|t}) = Var(\varepsilon_{t+1}) = \sigma^2
\tag{20}
$$

### 2기 앞 예측

2기 앞의 실제 데이터 생성 식은 (21)과 같다.

$$
Y_{t+2} = c + \phi Y_{t+1} + \varepsilon_{t+2} + \theta \varepsilon_{t+1}
\tag{21}
$$

#### 예측치 (Ŷ(t+2|t))

미래의 충격 $$\varepsilon_{t+2}, \varepsilon_{t+1}$$은 0이 되고, 미래 변수 $$Y_{t+1}$$은 앞서 구한 식 (18)의 1기 앞 예측치 $$\hat{Y}_{t+1|t}$$로 치환된다.

$$
\hat{Y}_{t+2|t} = c + \phi \hat{Y}_{t+1|t}
\tag{22}
$$

#### 예측오차 (e(t+2|t))

실제값에서 예측치를 뺀 후 식 (19)을 대입하여 정리한다.

$$
\begin{aligned}
e_{t+2|t} &= Y_{t+2} - \hat{Y}_{t+2|t} = \phi(Y_{t+1} - \hat{Y}_{t+1|t}) + \varepsilon_{t+2} + \theta \varepsilon_{t+1}  \\
&= \phi e_{t+1|t} + \varepsilon_{t+2} + \theta \varepsilon_{t+1}  \\
&= \phi \varepsilon_{t+1} + \varepsilon_{t+2} + \theta \varepsilon_{t+1}  \\
&= \varepsilon_{t+2} + (\phi + \theta)\varepsilon_{t+1} \qquad \text{(23)}
\end{aligned}
$$

#### 예측오차의 분산

&#8203;$$\varepsilon_{t+2}$$와 $$\varepsilon_{t+1}$$은 서로 독립이므로 교차항의 기댓값은 0이다.
따라서 각각의 분산을 더해준다.

$$
\begin{aligned}
Var(e_{t+2|t}) &= Var(\varepsilon_{t+2}) + (\phi + \theta)^2 Var(\varepsilon_{t+1})  \\
&= \sigma^2 + (\phi + \theta)^2 \sigma^2 = \sigma^2 \big[1 + (\phi + \theta)^2 \big] \qquad \text{(24)}
\end{aligned}
$$

### 3기 앞 예측

3기 앞의 실제 데이터 생성 식은 다음과 같다.

$$
Y_{t+3} = c + \phi Y_{t+2} + \varepsilon_{t+3} + \theta \varepsilon_{t+2}
\tag{25}
$$

#### 예측치 (Ŷ(t+3|t))

식 (25)에 조건부 기댓값을 취하고 식 (22)을 적용한다.

$$
\hat{Y}_{t+3|t} = c + \phi \hat{Y}_{t+2|t}
\tag{26}
$$

#### 예측오차 (e(t+3|t))

실제값에서 예측치를 빼어 1차적으로 정리한다.

$$
\begin{aligned}
e_{t+3|t} &= \phi(Y_{t+2} - \hat{Y}_{t+2|t}) + \varepsilon_{t+3} + \theta \varepsilon_{t+2}  \\
&= \phi e_{t+2|t} + \varepsilon_{t+3} + \theta \varepsilon_{t+2} \qquad \text{(27)}
\end{aligned}
$$

식 (27)에 앞서 구한 식 (23)의 $$e_{t+2|t}$$를 대입한다.

$$
\begin{aligned}
e_{t+3|t} &= \phi \big[ \varepsilon_{t+2} + (\phi + \theta)\varepsilon_{t+1} \big] + \varepsilon_{t+3} + \theta \varepsilon_{t+2}  \\
&= \varepsilon_{t+3} + (\phi + \theta)\varepsilon_{t+2} + \phi(\phi + \theta)\varepsilon_{t+1} \qquad \text{(28)}
\end{aligned}
$$

#### 예측오차의 분산

3개의 독립적인 오차항의 선형결합이므로, 각각의 계수를 제곱하여 더해준다.

$$
\begin{aligned}
Var(e_{t+3|t}) &= Var(\varepsilon_{t+3}) + (\phi + \theta)^2 Var(\varepsilon_{t+2}) + \phi^2(\phi + \theta)^2 Var(\varepsilon_{t+1})  \\
&= \sigma^2 + (\phi + \theta)^2 \sigma^2 + \phi^2(\phi + \theta)^2 \sigma^2  \\
&= \sigma^2 \big[ 1 + (\phi + \theta)^2 + \phi^2(\phi + \theta)^2 \big] \qquad \text{(29)}
\end{aligned}
$$

Horizon이 길어질수록 누적되는 오차의 분산이 점차 증가함을 식 (20), (24), (29)를 통해 확인할 수 있다.

## ARMA(p,q)

왈드 모형을 분석하려면 무한한 매개변수를 추정해야 하는데, 이는 불가능하다. 이를 위해, 계수 $$\psi_j$$들이 완전히 무작위로 분포하는 것이 아니라 특정한 규칙을 가지고 지수적으로 감소한다고 가정해보자.
가장 단순한 규칙으로 $$\psi_j = \phi^j$$ (단, $$|\phi| < 1$$) 라고 가정하자. Lag operator $$L$$을 사용해 표현하면 $$Y_t = (1 + \phi L + \phi^2 L^2 + \phi^3 L^3 + \cdots) \varepsilon_t$$가 된다.
&#8203;$$|\phi| < 1$$이므로, 괄호 안 무한등비급수는 수렴하여 $$\frac{1}{1 - \phi L}$$ 로 축약된다.

$$
Y_t = \frac{1}{1 - \phi L} \varepsilon_t
$$

위 식의 양변에 $$(1 - \phi L)$$을 곱하면 무한 역변환(invertibility)이 이루어져 다음 식이 도출된다.

$$
(1 - \phi L) Y_t = \varepsilon_t \implies Y_t = \phi Y_{t-1} + \varepsilon_t
$$

### p,q 결정 기준

#### AIC

$$
AIC = -2 \ln(\hat{L}) + 2k
$$

&#8203;$$\hat{L}$$은 Maximum Likelihood의 최댓값이며, $$k$$는 추정된 매개변수의 총 개수($$p+q+1$$ 등)이다.
만약 오차항이 정규분포를 따른다고 가정하고 RSS을 이용하여 표현하면 다음과 같이 근사할 수 있다 ($$T$$는 관측치 수).

$$
AIC = T \ln\left(\frac{RSS}{T}\right) + 2k
$$

- AIC는 새로운 데이터에 대한 예측 오차를 최소화하는 데 목적이 있다.
- 표본 크기 $$T$$가 무한히 커지더라도 true model을 선택할 확률이 1로 수렴하지 않는 inconsistency을 가진다. 즉, 표본이 커져도 불필요한 매개변수를 포함한 과대적합된 모형을 선택하는 경향이 있다.

#### SBC

슈바르츠가 posterior probability을 극대화하는 관점에서 도출한 기준으로, BIC (Bayesian Information Criterion)라고도 불린다.

$$
SBC = -2 \ln(\hat{L}) + k \ln(T)
$$

RSS을 이용한 형태는 다음과 같다.

$$
SBC = T \ln\left(\frac{RSS}{T}\right) + k \ln(T)
$$

- 매개변수 추가에 따른 페널티 항이 $$2k$$인 AIC와 달리, SBC는 $$k \ln(T)$$이다. 보통 관측치 수 $$T$$가 8 이상이면 $$\ln(T) > 2$$가 되므로, SBC가 AIC보다 모형의 복잡도에 훨씬 더 무거운 페널티를 부여한다.
- 표본 크기 $$T$$가 무한히 커질 때 참인 모형을 찾아낼 확률이 1로 수렴하는 일치성을 갖는다.

## Lag Operator

시차연산자 $$L$$은 임의의 시계열 변수 $$y_t$$의 시간 인덱스를 1기 후퇴시키는 linear operator로 정의된다.

$$
L y_t = y_{t-1}, \quad L^k y_t = y_{t-k} \quad (k \in \mathbb{N})
$$

시차다항식 $$\alpha(L) = 1 - \alpha L$$에 대하여, $$|\alpha| < 1$$을 만족할 때 이 연산자의 역연산자(inverse operator)는 무한 등비급수 형태로 전개된다.

$$
(1 - \alpha L)^{-1} = \sum_{i=0}^{\infty} \alpha^i L^i = 1 + \alpha L + \alpha^2 L^2 + \dots
$$

### AR(p) 모형의 정상성

&#8203;$$p$$차 자기회귀 모형 $$y_t = \phi_1 y_{t-1} + \dots + \phi_p y_{t-p} + \varepsilon_t$$는 시차연산자를 통해 다음과 같이 표현된다.

$$
\Phi(L) y_t = \varepsilon_t, \quad \text{where } \Phi(L) = 1 - \phi_1 L - \phi_2 L^2 - \dots - \phi_p L^p
$$

이 모형이 공분산 정상성을 만족하여 $$y_t = \Phi(L)^{-1}\varepsilon_t$$의 $$\text{MA}(\infty)$$ 형태로 전환되기 위한 필요충분조건은, 복소수 $$z$$에 대한 특성방정식 $$\Phi(z) = 0$$의 모든 근 $$z_i$$ ($$i=1, \dots, p$$)가 복소평면의 단위원 외부에 존재하는 것이다.

$$
\Phi(z) = 1 - \phi_1 z - \dots - \phi_p z^p = 0 \implies |z_i| > 1 \quad \forall i
$$

### MA(q) 모형의 가역성

&#8203;$$q$$차 이동평균 모형 $$y_t = \varepsilon_t + \theta_1 \varepsilon_{t-1} + \dots + \theta_q \varepsilon_{t-q}$$는 유한한 백색잡음의 선형결합이므로 파라미터에 관계없이 항상 정상적이다. 이 모형을 과거 관측치들의 무한합인 $$\text{AR}(\infty)$$ 형태, 즉 $$\Theta(L)^{-1}y_t = \varepsilon_t$$로 표현하여 현재의 오차항을 관측치로 역산해낼 수 있는 ‘가역성’을 만족하려면 다음 조건이 요구된다.

$$
\Theta(z) = 1 + \theta_1 z + \dots + \theta_q z^q = 0 \implies |z_i| > 1 \quad \forall i
$$

### ARMA(p, q) 모형의 종합 판별 및 공통근 상쇄

&#8203;$$\Phi(L)y_t = \Theta(L)\varepsilon_t$$ 모형에서 다음의 세 가지 검증이 동시에 이루어져야 한다.

**① 정상성:** $$\Phi(z) = 0$$의 모든 특성근이 $$|z| > 1$$을 만족하는가?

**② 가역성:** $$\Theta(z) = 0$$의 모든 특성근이 $$|z| > 1$$을 만족하는가?

**③ 식별성(Identification):** 다항식 $$\Phi(L)$$과 $$\Theta(L)$$ 사이에 공통근(Common root)이 존재하는가?

만약 $$\Phi(L) = (1 - \rho L)\Phi^*(L)$$이고 $$\Theta(L) = (1 - \rho L)\Theta^*(L)$$이라면, 양변에서 $$(1-\rho L)$$이 약분되어 실제 모형의 차수는 $$\text{ARMA}(p-1, q-1)$$로 축소된다.

## 상태공간 표현과 동반행렬(Companion Matrix) 분석

고차 차분방정식인 $$\text{AR}(p)$$ 모형의 동학을 선형대수적으로 취급하기 위해, $$p \times 1$$ 상태 벡터 $$\boldsymbol{Y}_t$$와 $$p \times p$$ 동반행렬 $$\boldsymbol{F}$$를 도입하여 1계 벡터 시스템($$\text{VAR}(1)$$)으로 압축한다.

$$
\underbrace{
\begin{bmatrix}
y_t \\ y_{t-1} \\ y_{t-2} \\ \vdots \\ y_{t-p+1}
\end{bmatrix}
}_{\boldsymbol{Y}_t}
=
\underbrace{
\begin{bmatrix}
\phi_1 & \phi_2 & \dots & \phi_{p-1} & \phi_p \\
1 & 0 & \dots & 0 & 0 \\
0 & 1 & \dots & 0 & 0 \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & \dots & 1 & 0
\end{bmatrix}
}_{\boldsymbol{F}}
\underbrace{
\begin{bmatrix}
y_{t-1} \\ y_{t-2} \\ y_{t-3} \\ \vdots \\ y_{t-p}
\end{bmatrix}
}_{\boldsymbol{Y}_{t-1}}
+
\underbrace{
\begin{bmatrix}
\varepsilon_t \\ 0 \\ 0 \\ \vdots \\ 0
\end{bmatrix}
}_{\boldsymbol{\mathcal{E}}_t}
$$

즉, $$\boldsymbol{Y}_t = \boldsymbol{F}\boldsymbol{Y}_{t-1} + \boldsymbol{\mathcal{E}}_t$$가 성립한다.

#### Eigenvalue 기반의 안정성 조건

이 벡터 시스템이 발산하지 않고 안정적으로 수렴하기 위한 필요충분조건은 동반행렬 $$\boldsymbol{F}$$의 모든 고유값 $$\lambda_i$$ ($$i=1, \dots, p$$)의 모듈러스가 단위원 내부에 존재하는 것이다.

$$
\det(\lambda \boldsymbol{I} - \boldsymbol{F}) = 0 \implies |\lambda_i| < 1 \quad \forall i
$$

#### 특성근(z)과 고유값(λ)의 쌍대성 증명

시차다항식의 특성방정식 $$\Phi(z)=0$$과 동반행렬의 특성방정식 $$\det(\lambda \boldsymbol{I} - \boldsymbol{F}) = 0$$의 해는 서로 역수 관계($$\lambda = z^{-1}$$)에 놓인다.

동반행렬 $$\boldsymbol{F}$$의 특성방정식을 Laplace expansion으로 전개하면 다음과 같은 $$p$$차 다항식을 얻는다.

$$
\det(\lambda \boldsymbol{I} - \boldsymbol{F}) = \lambda^p - \phi_1 \lambda^{p-1} - \phi_2 \lambda^{p-2} - \dots - \phi_p = 0
$$

위 식의 $$\lambda$$가 0이 아님을 전제하고 양변을 $$\lambda^p$$로 나눈다.

$$
1 - \phi_1 \left(\frac{1}{\lambda}\right) - \phi_2 \left(\frac{1}{\lambda}\right)^2 - \dots - \phi_p \left(\frac{1}{\lambda}\right)^p = 0
$$

여기서 $$z = \lambda^{-1}$$로 치환하면, 정확히 시차연산자를 통해 도출했던 특성방정식 $$\Phi(z) = 0$$과 일치한다.

$$
1 - \phi_1 z - \phi_2 z^2 - \dots - \phi_p z^p = 0 \quad \iff \quad \Phi(z) = 0
$$

따라서 임의의 $$i$$번째 근에 대하여 $$\lambda_i = \frac{1}{z_i}$$이 성립하므로, 다음의 두 안정성 판별 기준은 동치이다.

$$
|z_i| > 1 \quad \text{(Outside Unit Circle)} \iff |\lambda_i| < 1 \quad \text{(Inside Unit Circle)}
$$
