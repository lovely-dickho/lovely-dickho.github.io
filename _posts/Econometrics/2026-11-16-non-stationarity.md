---
title: "[Econometrics] Non-Stationarity"
date: 2026-11-16 18:57:00 +0900
categories: [Econometrics]
tags: [비정상성, 단위근, ADF, 가성회귀, 랜덤워크]
math: true
toc: true
---

## 시간추세와 확률추세

#### 시간추세 (확정추세)

시간에 의해 결정되는 추세이다. 시간추세 모형의 단순한 형태로 다음을 들 수 있다.

$$
Z_t = \beta_0 + \beta_1 t + u_t
$$

단, $$t$$는 시간지수(time index)이고 $$u_t$$는 정상 시계열변수이다. $$\beta_0 + \beta_1 t$$는 시간지수 $$t$$에 의해 결정되는 시간추세를 나타낸다. $$u_t$$는 ARMA($$p,q$$) 또는 백색잡음과 같은 정상과정을 따르는 확률과정을 나타낸다.

#### 확률추세

일단 무작위적 충격이 발생하면 그것의 영향이 영속적으로 지속되기 때문에 확률적으로 형성된 추세를 의미한다.

### Random Walk Process

RW process는 가장 단순한 형태로서, 확률추세를 갖는 확률과정의 전형적 특성을 보여준다.

$$
Y_t = Y_{t-1} + \varepsilon_t \quad \varepsilon_t \sim N(0, \sigma_\varepsilon^2)
$$

#### 평균과 분산

RW은 그것의 분산이 시간 $$t$$에 대해 일정하지 않기 때문에 비정상과정이다. 0시점의 초기값을 $$y_0$$로 가정하고, 1시점에서 $$t$$시점까지 전방적 반복계산을 해보자.

$$
\begin{aligned}
&Y_1 = y_0 + \varepsilon_1 \\
&Y_2 = Y_1 + \varepsilon_2 = y_0 + \varepsilon_1 + \varepsilon_2 \\
&Y_3 = Y_2 + \varepsilon_3 = y_0 + \varepsilon_1 + \varepsilon_2 + \varepsilon_3
\end{aligned}
$$

그 결과, 다음을 얻는다.

$$
Y_t = y_0 + \varepsilon_1 + \varepsilon_2 + \cdots + \varepsilon_t
$$

이 식의 우변에서 모든 $$\varepsilon$$이 0의 평균을 갖는 백색잡음이므로, 기대치 $$\mathbb{E}(Y_t)$$는 모든 $$t$$에 대하여 $$y_0$$로서 일정한 값이다. $$Y_t$$의 분산은 다음과 같이 구해진다.

$$
\sigma_t^2 = \mathbb{E}(\varepsilon_1 + \varepsilon_2 + \cdots + \varepsilon_t)^2 = t\sigma_\varepsilon^2
$$

분산은 시간종속적이며 $$t$$가 증가하면 분산도 증가한다. 무작위적 충격이 시계열과정에 미치는 영향이 영속적으로 지속되고, 그러한 무작위적 충격이 누적됨에 따라 $$t$$가 증가할수록 분산도 증가한다. 따라서 RW process는 비정상과정이며, 정상적 과정에만 사용되는 ARMA(p,q)모형으로 모형화될 수 없다.

#### k시차 자기상관성

$$
\begin{aligned}
\rho_k &= \frac{\mathbb{E}[(\varepsilon_t + \cdots + \varepsilon_1)(\varepsilon_{t-k} + \cdots + \varepsilon_1)]}{\sqrt{t\sigma_\varepsilon^2}\sqrt{(t-k)\sigma_\varepsilon^2}} \\
&= \frac{(t-k)\sigma_\varepsilon^2}{\sqrt{t\sigma_\varepsilon^2}\sqrt{(t-k)\sigma_\varepsilon^2}} = \frac{\sqrt{t-k}}{\sqrt{t}}
\end{aligned}
$$

만일 $$k$$가 $$t$$에 비해 작은 값이라면, AC는 1에 가까워진다.

#### 예측

&#8203;$$\{y_1, \cdots, y_T\}$$의 과거자료가 주어져 있을 때, RW process의 $$Y_{T+l}$$ 값을 얻기 위해 전방적 반복계산을 행하면 다음 식이 얻어진다.

$$
Y_{T+l} = y_T + \varepsilon_{T+1} + \varepsilon_{T+2} + \cdots + \varepsilon_{T+l}
$$

RW process의 $$l$$기후 예측치는 $$Y_{T+l}$$의 조건부기대치이므로, $$l$$기후 예측치는 다음과 같다.

$$
\mathbb{E}_T(Y_{T+l}) = y_T
$$

즉, RW process의 $$l$$기후 예측치는 예측기간에 관계없이 예측원점인 $$T$$시점의 관찰치와 같다.

#### 시차연산자와 RW

lag operator를 사용하여 RW과정을 나타내면 다음과 같다.

$$
Y_t - Y_{t-1} = (1 - L)Y_t = \varepsilon_t
$$

이 식은 RW의 1차 차분이 백색잡음과정임을 보여준다. 비정상과정은 차분함으로써 정상과정을 얻을 수 있다. 이와 같이 한 번만 차분하여 정상과정이 얻어지는 시계열과정은 **‘1차 적분되었다(integrated of order 1)’**고 한다.

### 랜덤워크 플러스 드리프트 과정(random walk plus drift)

$$
Y_t = \delta + Y_{t-1} + \varepsilon_t
$$

&#8203;$$\delta$$는 드리프트(drift)라고 부르는 상수항이다. 0시점의 초기값을 $$y_0$$라고 가정하면, RW 플러스 드리프트과정은 후방적 반복계산에 의해 다음 식으로 다시 쓸 수 있다.

$$
Y_t = y_0 + t\delta + \varepsilon_1 + \cdots + \varepsilon_t
$$

위 식은 RW 플러스 드리프트과정이 확정추세와 확률추세 모두를 갖고 있다는 것을 보여주고 있다. 두 번째 항은 시계열과정이 확정적으로 $$\delta$$의 기울기를 갖고 증가한다는 것을 의미한다. 이와 함께 시계열과정에 주어진 모든 무작위적 충격의 영향이 영속적으로 사라지지 않기 때문에 이 시계열과정은 확률추세를 갖는다.

### ARIMA(p,d,q)

어떤 시계열이 $$d$$차 차분이 필요되면 그 시계열은 $$d$$차 적분되었다고 말하며 ARIMA($$p,d,q$$)로 나타낸다.

ARIMA(1,1,1)과정은 1차 차분에 의해 ARMA(1,1)로 변환시킬 수 있는 시계열으로서 다음과 같이 표현된다.

$$
\underbrace{(1 - \phi L)}_{\text{AR}} \underbrace{(1 - L)}_{\text{MA}} Y_t = (1 + \theta L) \varepsilon_t
$$

위의 식에서 $$(1 - L)$$은 1차 적분되어 있다는 것을 의미한다. AR 부분의 시차 다항식인 $$(1 - \phi L)(1 - L)$$은 $$(1 - \phi L)$$과 $$(1 - L)$$로 인수분해된다. 즉, AR 부분의 시차 다항식의 해가 1이라는 것, 즉 단위근(unit root)을 갖고 있다는 것을 보여준다.

동일한 규칙에 의해, ARIMA($$p,d,q$$)과정은 다음과 같다.

$$
(1 - \phi_1 L - \cdots - \phi_p L^p)(1 - L)^d Y_t = (1 + \theta_1 L + \cdots + \theta_q L^q) \varepsilon_t
$$

시계열과정의 AR 부분이 단위근을 갖고 있다. 이처럼 단위근의 존재는 확률추세를 갖고 있는 비정상 과정의 매우 중요한 특성이다. 예를 들어, RW process는 ARIMA($$p,d,q$$)에서 $$d=1, p=q=0$$인 경우로서 단위근을 갖고 있다.

## 가성회귀와 공적분

#### 가성회귀(Spurious Regression)와 그 원인

가성회귀는 아무런 인과관계를 갖지 않는 변수들의 회귀분석이 높은 적합도를 갖는 경우를 의미한다. 이때의 회귀관계는 높은 적합도에도 불구하고 어떤 진정한 경제적 의미를 갖고 있지 않다.

- 가성회귀는 Granger & Newbold(1974)가 몬테카를로 시뮬레이션을 사용하여 처음 발견한 것으로서, 서로 독립적이고 인과관계를 갖지 않는 둘 또는 그 이상의 비정상적(nonstationary) 변수들에 의해 만들어진 회귀관계를 말한다.
- 이는 오차항 $$u_t$$가 RW process와 동일한 특성을 갖고 있기 때문이다. $$u_t$$에 주어진 무작위적 충격이 소멸하지 않고, $$t$$가 증가함에 따라 $$u_t$$의 분산도 증가한다. $$u_t$$는 실질적으로 1시차에서 높은 자기상관성을 갖고 있는 RW process이다.

#### 가성회귀의 처방

가성회귀의 문제에 대한 처방으로는 종속변수와 독립변수의 자료를 차분(differencing)하고, 차분된 자료를 사용하여 OLS를 적용하거나 Cochrane-Orcutt 방법을 사용하는 것이 있다.

#### Cointegration

가성회귀의 한 가지 특별한 경우는 공적분된 경우이다. **만일 비정상 시계열과정 $$X_t$$와 $$Y_t$$가 동일 차수로 적분되어 있고, 두 시계열과정의 선형결합이 정상과정을 만든다면, 우리는 $$X_t$$와 $$Y_t$$가 '공적분되었다(cointegrated)' 라고 말한다.**

### Cochrane-Orcutt 방법

코크란-오르컷 방법은 선형회귀모형의 오차항에 1차 자기상관이 존재할 때 이를 교정하여 효율적인 추정량을 얻기 위한 반복적인 FGLS의 일종이다.

고전적 선형회귀모형에서 오차항이 1차 자기상관, 즉 AR(1) 과정을 따른다고 가정하자.

$$
\begin{aligned}
Y_t &= \beta_0 + \beta_1 X_t + u_t \\
u_t &= \rho u_{t-1} + \varepsilon_t, \quad (|\rho| < 1)
\end{aligned}
$$

여기서 $$\varepsilon_t$$는 고전적 가정을 충족하는 백색잡음이다. 위 모형의 양변에 $$\rho$$를 곱하고 1기 지연시키면 다음과 같다.

$$
\rho Y_{t-1} = \rho \beta_0 + \rho \beta_1 X_{t-1} + \rho u_{t-1}
$$

원래의 식에서 지연된 식을 차감하면 다음과 같은 변환된 회귀식을 얻을 수 있다.

$$
\underbrace{Y_t - \rho Y_{t-1}}_{Y_t^*} = \underbrace{\beta_0(1 - \rho)}_{\beta_0^*} + \beta_1 \underbrace{(X_t - \rho X_{t-1})}_{X_t^*} + \underbrace{(u_t - \rho u_{t-1})}_{\varepsilon_t}
$$

&#8203;$$\varepsilon_t$$는 자기상관이 제거된 상태이므로 OLS를 적용하여 통계적으로 유의미한 추정량을 얻을 수 있다.

실제로는 자기상관계수 $$\rho$$를 알 수 없으므로, 다음의 과정을 통해 모형을 반복적으로 추정한다.

**①** 원래의 모형 $$Y_t = \beta_0 + \beta_1 X_t + u_t$$를 OLS로 추정하여 잔차 $$\hat{u}_t$$를 구한다.

**②** 구해진 잔차를 이용하여 $$\hat{u}_t = \rho \hat{u}_{t-1} + v_t$$ 형태의 보조회귀를 수행함으로써 $$\hat{\rho}$$를 추정한다.

**③** 추정된 $$\hat{\rho}$$를 이용하여 변수들을 $$Y_t^* = Y_t - \hat{\rho} Y_{t-1}$$ 및 $$X_t^* = X_t - \hat{\rho} X_{t-1}$$로 준차분(quasi-differencing) 변환한다.

**④** 변환된 종속변수 $$Y_t^*$$를 변환된 독립변수 $$X_t^*$$에 대해 회귀분석하여 새로운 모수 추정치 $$\hat{\beta}_0, \hat{\beta}_1$$을 구한다.

**⑤** 새로운 모수 추정치를 바탕으로 원래 모형에서 다시 새로운 잔차를 구하고, $$\hat{\rho}$$를 재추정한다. 이 과정을 $$\hat{\rho}$$ 값이나 추정된 모수들이 수렴할 때까지(이전 단계의 추정치와 유의미한 차이가 없을 때까지) 반복한다.

## 추세의 제거

추세에는 시간추세(time trend)와 확률추세(stochastic trend)의 두 종류가 있다. 이 두 가지 추세를 제거할 때에는 각기 다른 방법을 사용하여야 한다.

### Detrending

시간추세는 시간에 의해 만들어지는 것이므로, 이를 제거하기 위해서는 디트렌딩해야 한다. 디트렌딩은 시계열 자료를 시간에 대하여 회귀분석함으로써 시간추세가 제거된 잔차 시계열을 구하는 방법이다.

만일 관찰된 표본 $$\{y_1, \cdots, y_T\}$$가 선형 시간추세를 갖고 있다고 생각되면, 다음의 모형을 회귀분석할 수 있다.

$$
y_t = \beta_0 + \beta_1 t + u_t
$$

시간추세가 제거된 잔차 $$\hat{u}_t$$는 $$\hat{u}_t = y_t - (\hat{\beta}_0 + \hat{\beta}_1 t)$$에 의해 계산된다. 시간추세가 제거된 잔차 $$\hat{u}_t$$가 예상대로 정상 시계열이고 ARMA($$p,q$$)과정을 갖는다면, 원래의 시계열인 $$Y_t$$는 추세정상적(trend-stationary : TS) 과정이라고 한다.

### Differencing

확률추세는 그 영향이 영속적으로 지속되는 충격의 누적에 의해 만들어진다. 확률추세를 제거하기 위해서는 차분이 사용된다. 대부분의 비정상 시계열들은 1차 차분만을 필요로 하지만, 예외적인 경우 고차 차분이 필요한 경우도 있다.

확률추세를 가진 시계열을 차분함으로써 정상과정이 얻어질 수 있다면, 우리는 그런 시계열을 차분정상적(difference-stationary : DS) 과정이라고 부른다. 예를 들어, $$Y_t$$를 차분하여 정상 시계열인 ARMA(1,1) 시계열이 만들어진다면, $$Y_t$$는 DS 시계열이다. 이런 시계열과정의 AR 부분의 다항식은 단위근을 갖기 때문에, 차분정상적(DS) 과정은 단위근과정(unit root process)이라고도 부른다.

### 디트렌딩 vs 차분

정상과정을 얻기 위해 TS 과정을 디트렌딩하지 않고 차분하면 어떤 문제가 발생하는지 보자.

$$
Z_t = \beta_0 + \beta_1 t + \varepsilon_t \quad \varepsilon_t \sim N(0, \sigma_\varepsilon^2)
$$

&#8203;$$Z_t$$의 1차 차분은 다음과 같다.

$$
\begin{aligned}
(1 - L)Z_t &= Z_t - Z_{t-1} = \beta_0 + \beta_1 t + \varepsilon_t - [\beta_0 + \beta_1 (t-1) + \varepsilon_{t-1}] \\
&= \beta_1 + \varepsilon_t - \varepsilon_{t-1}
\end{aligned}
$$

MA 부분은 $$(1 - L)\varepsilon_t$$로서 단위근을 갖고 있으므로, $$Z_t$$의 1차 차분은 가역적(invertible)이지 않다. 따라서 차분된 시계열을 사용하여 정상적 모형을 추정할 수 없다.

DS 과정을 차분하지 않고 디트렌딩해보자.

$$
\begin{aligned}
Y_t &= \delta + Y_{t-1} + \varepsilon_t \\
&= y_0 + t\delta + \varepsilon_1 + \cdots + \varepsilon_t
\end{aligned}
$$

RW 플러스 드리프트 모형은 $$t\delta$$항 때문에 TS과정과 매우 유사하게 보인다. 만일 디트렌딩을 하기 위해 $$t$$에 대해 $$Y_t$$를 회귀분석하고 $$\delta$$를 추정하면, 시간추세를 제거한 시계열이 다음과 같이 얻어진다.

$$
Y_t - t\delta = y_0 + \varepsilon_1 + \cdots + \varepsilon_t
$$

시간추세를 제거한 시계열은 RW process와 정확하게 동일하다. 그러므로 RW 플러스 드리프트과정을 디트렌딩하게 되면 정상 과정이 얻어지지 않는다.

## 단위근 검증

### ADF 검증

DS 과정의 가장 간단한 형태는 RW process이고, RW process는 AR(1)과정 중 단위근을 갖고 있는 특별한 경우이다. AR(1)모형, $$Y_t = \phi Y_{t-1} + \varepsilon_t$$의 두 변에서 $$Y_{t-1}$$을 빼면 다음 식이 얻어진다.

$$
\Delta Y_t = (\phi - 1) Y_{t-1} + \varepsilon_t = \gamma Y_{t-1} + \varepsilon_t
$$

단위근검증은 $$H_0 : \gamma = 0$$을 검증하는 것으로서, $$H_0 : \phi = 1$$의 검증과 동일한 것이다. 대립가설은 $$H_a : \phi < 1$$ 또는 $$H_a : \gamma < 0$$이다. 단위근검증은 단측 가설검증으로, $$\gamma$$의 추정치가 0과 같거나 크면, 그 시계열이 단위근을 갖는다고 결론내린다.

#### Augmented Dickey-Fuller(ADF) tests

대부분의 시계열은 drift, 시간추세 또는 AR 부분을 갖고 있다. 이를 모형화하면 다음과 같다.

$$
\begin{aligned}
\Delta Y_t &= \gamma Y_{t-1} + \alpha_1 \Delta Y_{t-1} + \cdots + \alpha_p \Delta Y_{t-p} + \varepsilon_t \qquad \text{(1)}\\
\Delta Y_t &= \delta + \gamma Y_{t-1} + \alpha_1 \Delta Y_{t-1} + \cdots + \alpha_p \Delta Y_{t-p} + \varepsilon_t \qquad \text{(2)} \\
\Delta Y_t &= \delta + \beta t + \gamma Y_{t-1} + \alpha_1 \Delta Y_{t-1} + \cdots + \alpha_p \Delta Y_{t-p} + \varepsilon_t \qquad \text{(3)}
\end{aligned}
$$

여기서 $$\varepsilon_t$$는 백색잡음이다. 이 세 가지 모형에는 차분항의 시차값인 $$\alpha_i \Delta Y_{t-i}$$가 덧붙여져 있다.

- $$Y_t$$가 드리프트와 시간추세를 갖지 않다고 생각되는 적분과정이면, (1)를 사용
- 만일 $$Y_t$$가 오직 드리프트만을 갖고 있다고 생각되면 (2)를 사용
- $$Y_t$$가 드리프트와 시간추세를 갖고 있다고 판단될 때에는, (3)를 사용

**중요: 검증 대상은 $$\gamma$$이며, 세 가지 경우 모두에 대하여 검증하고자 하는 귀무가설은 $$H_0 : \gamma = 0$$이다.**

#### 주의

Dickey & Fuller는 표본크기가 동일하다 하더라도, 위의 세 가지 모형 중 어떤 것을 사용하느냐에 따라 ADF 검증통계량의 임계값이 크게 달라진다는 것을 발견하였다.

#### 단위근 검증 전략

**①** 절편항과 시간추세를 갖는 (3)를 추정한다. 그리고 다음 단계로 나가기 위해 다음 세 가지 중 하나를 선택한다.

- 만일 $$\gamma$$가 0이 아니라면 시계열에 단위근이 없다고 결론내린다.
- 만일 $$\gamma$$와 $$\beta$$ 모두 0이라면 2로 진행한다.
- 만일 $$\gamma$$가 0이고 $$\beta$$가 0이 아니라면, 주어진 시계열은 시간추세와 단위근을 갖고 있다고 결론내린다.

**②** 드리프트가 있는 DS과정을 가정하고 있는 (2)를 추정한다. 여기에서도 다음 세 가지 중 하나를 선택한다.

- 만일 $$\gamma$$가 0이 아니라면 단위근이 없다고 결론내린다.
- 만일 $$\gamma$$와 $$\delta$$ 모두 0이 아니라면 3으로 진행한다.
- 만일 $$\gamma$$가 0이고 $$\delta$$가 0이라면, 주어진 시계열은 드리프트를 갖고 있고 동시에 단위근을 갖고 있다고 결론내린다.

**③** 순수한 DS과정을 가정하고 있는 (1)를 추정한다. 여기에서는 다음 둘 중 하나를 선택한다.

- 만일 $$\gamma$$가 0이 아니라면 단위근이 없다고 결론내린다.
- 만일 $$\gamma$$가 0이라면 주어진 시계열이 단위근을 갖고 있다고 결론내린다.
