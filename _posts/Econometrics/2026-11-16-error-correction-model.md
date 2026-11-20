---
title: "[Econometrics] Error Correction Model"
date: 2026-11-16 18:56:00 +0900
categories: [Econometrics]
tags: [ECM, VECM, 공적분, Granger]
math: true
toc: true
---

## ECM

I(1) 변수들 사이의 관계를 장기관계, I(0) 변수들 사이의 관계를 단기관계라 한다. **단기 오차 수정 모형**은 공적분 관계를 끼워 넣은 I(0) 변수들 사이의 동태적 관계를 다룬다.

ARDL(1, 1) 모형 $$y_t = \delta + \theta_1 y_{t-1} + \delta_0 x_t + \delta_1 x_{t-1} + v_t$$에 대한 ECM은 다음과 같다.

$$
\begin{aligned}
\Delta y_t &= -\alpha (y_{t-1} - \beta_1 - \beta_2 x_{t-1}) + \delta_0 \Delta x_t + v_t
\end{aligned}
$$

위 식은 다음과 같은 이유로 **오차 수정식**이라고 한다.

- $$y_{t-1} - \beta_1 - \beta_2 x_{t-1}$$은 $$y_{t-1}$$의 장기적인 값, $$\beta_1 + \beta_2 x_{t-1}$$로부터 $$y_{t-1}$$이 벗어난 정도, 다시 말해 이전 기간에서의 ‘오차’를 보여준다.
- $$(\theta_1 - 1)$$항은 ‘오차’에 대한 $$\Delta y_t$$의 ‘수정’을 보여준다. $$y$$와 $$x$$ 사이에 공적분 관계가 존재해서 ‘오차가 수정되도록’ 조정이 언제나 작동하는 경우, 경험적으로 볼 때 $$\theta_1 < 1$$이라는 사실을 의미하는 $$(1 - \theta_1) > 0$$을 또한 발견해야만 한다.

### ECM 도출

ARDL(1, 1) 모형 $$y_t = \delta + \theta_1 y_{t-1} + \delta_0 x_t + \delta_1 x_{t-1} + v_t$$에서 장기 정상상태를 가정하면 장기균형관계를 얻을 수 있다.

$$
y^* = \beta_1 + \beta_2 x^* \quad where \quad \beta_1 = \frac{\delta}{1 - \theta_1}, \beta_2 = \frac{\delta_0 + \delta_1}{1 - \theta_1}
$$

이를 바탕으로 오차수정모형(ECM)을 도출하는 두 가지 방법을 정리하면 다음과 같다.

#### 원래 방식

기본 ARDL 모형의 양변에 $$y_{t-1}$$을 빼고, 우변에 $$\delta_0 x_{t-1}$$을 더하고 빼서 식을 변형한다.

$$
\begin{aligned}
\Delta y_t &= \delta + (\theta_1 - 1) y_{t-1} + \delta_0 x_t + \delta_1 x_{t-1} + v_t \\
&= \delta + (\theta_1 - 1) y_{t-1} + \delta_0 \Delta x_t + (\delta_0 + \delta_1) x_{t-1} + v_t
\end{aligned}
$$

우변의 $$y_{t-1}$$ 계수인 $$(\theta_1 - 1)$$로 관련 항들을 묶어낸다.

$$
\begin{aligned}
\Delta y_t &= (\theta_1 - 1) \left[ \frac{\delta}{\theta_1 - 1} + y_{t-1} + \frac{\delta_0 + \delta_1}{\theta_1 - 1} x_{t-1} \right] + \delta_0 \Delta x_t + v_t
\end{aligned}
$$

괄호 안의 식을 정상상태 모수 $$\beta_1, \beta_2$$ 형태로 변환하고, 조정속도 파라미터 $$\alpha = 1 - \theta_1$$을 적용하면 최종 ECM이 도출된다.

$$
\begin{aligned}
\Delta y_t &= -\alpha (y_{t-1} - \beta_1 - \beta_2 x_{t-1}) + \delta_0 \Delta x_t + v_t
\end{aligned}
$$

#### 오차 정의

장기균형 오차 $$e_{t-1} = y_{t-1} - \beta_1 - \beta_2 x_{t-1}$$를 먼저 정의한다. 이 식의 양변에 $$-\alpha$$ (즉, $$\theta_1 - 1$$)를 곱하여 원래의 ARDL 모수 형태로 전개한다.

$$
\begin{aligned}
-\alpha e_{t-1} &= (\theta_1 - 1) \left( y_{t-1} - \frac{\delta}{1 - \theta_1} - \frac{\delta_0 + \delta_1}{1 - \theta_1} x_{t-1} \right) \\
-\alpha e_{t-1} &= (\theta_1 - 1)y_{t-1} + \delta + (\delta_0 + \delta_1)x_{t-1}
\end{aligned}
$$

이를 $$\delta + (\theta_1 - 1)y_{t-1}$$ 에 대해 정리한다.

$$
\begin{aligned}
\delta + (\theta_1 - 1)y_{t-1} &= -\alpha e_{t-1} - \delta_0 x_{t-1} - \delta_1 x_{t-1}
\end{aligned}
$$

양변에서 $$y_{t-1}$$을 뺀 ARDL 모형($$\Delta y_t = \delta + (\theta_1 - 1) y_{t-1} + \delta_0 x_t + \delta_1 x_{t-1} + v_t$$)에 위 식을 통째로 대입하여 정리한다.

$$
\begin{aligned}
\Delta y_t &= \left[ -\alpha e_{t-1} - \delta_0 x_{t-1} - \delta_1 x_{t-1} \right] + \delta_0 x_t + \delta_1 x_{t-1} + v_t \\
\Delta y_t &= -\alpha e_{t-1} - \delta_0 x_{t-1} + \delta_0 x_t + v_t \\
\Delta y_t &= -\alpha e_{t-1} + \delta_0 \Delta x_t + v_t
\end{aligned}
$$

## 공적분

공적분 분석은 **공적분 벡터의 존재와 수를 검증하는 것**과 **공적분 벡터를 추정하는 것**으로 이루어져 있다.

### 공적분 벡터

1차적분되는 세 변수 $$X_{1,t}, X_{2,t}, X_{3,t}$$가 있다고 하자. 이 변수들의 선형결합은 다음과 같다.

$$
e_t = \beta_1 X_{1,t} + \beta_2 X_{2,t} + \beta_3 X_{3,t}
\tag{1}
$$

- 선형결합 $$e_t$$가 정상과정을 따른다면, 변수 $$X_{1,t}, X_{2,t}, X_{3,t}$$는 (1,1)차 공적분되어 있다고 말하고 $$CI(1,1)$$로 표시한다. (1)은 공적분방정식이라고 한다.
- 선형결합의 계수벡터 $$\boldsymbol{\beta} = (\beta_1, \beta_2, \beta_3)^T$$는 공적분벡터라 한다. 공적분벡터는 정규화될 수 있다.

  $$
  \boldsymbol{\beta}^* = \left(1, \frac{\beta_2}{\beta_1}, \frac{\beta_3}{\beta_1}\right)^T
  $$

  $$\boldsymbol{\beta}^*$$ 는 선형의존적 벡터 $$\boldsymbol{\beta}$$의 배수이다.

공적분벡터를 정규화하면, (1)의 선형결합은 다음과 같다.

$$
X_{1,t} = - \frac{\beta_2}{\beta_1} X_{2,t} - \frac{\beta_3}{\beta_1} X_{3,t} + \frac{1}{\beta_1} e_t
\tag{2}
$$

- (1)에서 선형결합 $$e_t$$는 (2)의 회귀모형에서 잔차항으로 나타나게 되고, 공적분벡터는 회귀계수들의 벡터로 나타나고 있다.
- (2)과 일반적인 회귀모형의 차이점은 (2)의 종속변수와 독립변수가 동일한 차수로 적분되어 있다는 것이다.
- 잔차항 $$e_t$$는 정상과정이다.

#### 경제적 의미

- 어떤 경제이론은 이 적분된 변수들 사이의 장기균형관계를 규정한다. 여기서 ‘균형’은 엄격한 의미의 최적 상태가 아니라 넓은 의미에서의 안정된 상태이다. 만일 그 이론들이 실제로 성립한다면, 변수들의 선형결합은 정상시계열과정을 가져야만 한다.
- (1)에서 변수들의 값이 정확한 장기균형을 이룬다면, $$e_t$$는 0이다. 그러나 예상치 못한 충격으로 인해 장기균형에서 일시적 또는 단기적 이탈이 발생하면, $$e_t$$는 0이 아니다. $$e_t$$가 정상적이라는 사실은 단기적 이탈이 오래 지속되지 않는다는 것을 의미한다.

#### 일반화

어떤 변수벡터 $$\boldsymbol{X} = (X_{1,t}, \cdots, X_{k,t})^T$$의 성분들이 모두 $$d$$차 적분되어 있고 그것들을 선형결합한 것이 $$(d-b)$$차 적분된 시계열과정이라면, 변수 $$X_{1,t}, \cdots, X_{k,t}$$는 $$(d,b)$$차 공적분되어 있다고 하고 $$\boldsymbol{X} = CI(d,b)$$로 표시한다.

&#8203;$$k$$개의 적분된 변수들로 만들 수 있는 공적분벡터의 수는 최대 $$k-1$$개까지 가능하다.

- 모든 성분들이 1차 적분되는 3개의 변수들로 이루어진 변수벡터 $$\boldsymbol{X} = (X_{1,t}, X_{2,t}, X_{3,t})^T$$가 있다고 하자.
- 이 변수들이 다음의 공적분벡터를 가질 수 있다고 가정하자.

  $$
  \begin{aligned}
  \boldsymbol{\beta}_1 &= (\beta_{11}\ \beta_{21}\ \beta_{31})^T \\
  \boldsymbol{\beta}_2 &= (\beta_{12}\ \beta_{22}\ \beta_{32})^T
  \end{aligned}
  $$

- 변수의 선형결합은 다음과 같다.

  $$
  \begin{pmatrix}
  \beta_{11} & \beta_{21} & \beta_{31} \\
  \beta_{12} & \beta_{22} & \beta_{32}
  \end{pmatrix}
  \begin{pmatrix}
  X_{1,t} \\
  X_{2,t} \\
  X_{3,t}
  \end{pmatrix}
  =
  \begin{pmatrix}
  e_{1,t} \\
  e_{2,t}
  \end{pmatrix}
  \tag{3}
  $$

  선형결합 $$e_{1,t}$$와 $$e_{2,t}$$는 정상적이다.
- 공적분벡터로 구성된 행렬을 공적분행렬이라고 부른다. (3)에서 공적분행렬 $$\boldsymbol{\beta}$$는 다음과 같다.

  $$
  \boldsymbol{\beta} =
  \begin{pmatrix}
  \beta_{11} & \beta_{12} \\
  \beta_{21} & \beta_{22} \\
  \beta_{31} & \beta_{32}
  \end{pmatrix}
  $$

- 일반적으로, $$\boldsymbol{X} = (X_{1,t}, \cdots, X_{k,t})^T$$의 모든 성분들이 1차 적분되어 있고, $$r$$개(단, $$r \leq k-1$$)개의 공적분벡터가 존재한다면, 공적분행렬의 $$rank(\boldsymbol{\beta})$$는 $$r$$이다.

## VECM

&#8203;$$\boldsymbol{X} = (X_{1,t}, \cdots, X_{k,t})^T$$의 모든 성분들이 1차 적분되어 있다고 가정하자. $$\boldsymbol{X}$$는 다음의 VAR($$p$$)모형으로 표현될 수 있다.

$$
\boldsymbol{X}_t = \boldsymbol{A}_1 \boldsymbol{X}_{t-1} + \cdots + \boldsymbol{A}_p \boldsymbol{X}_{t-p} + \boldsymbol{\varepsilon}_t
\tag{4}
$$

단, $$\boldsymbol{X}_t$$는 ($$k \times 1$$) 변수벡터이고, $$\boldsymbol{A}_i$$는 ($$k \times k$$) 계수행렬이고, $$\boldsymbol{\varepsilon}_t$$는 백색잡음항들의 ($$k \times 1$$)벡터

이 VAR($$p$$)모형은 다음으로 변환될 수 있다.

$$
\Delta \boldsymbol{X}_t = \boldsymbol{\Pi} \boldsymbol{X}_{t-1} + \sum_{i=1}^{p-1} \boldsymbol{\Gamma}_i \Delta \boldsymbol{X}_{t-i} + \boldsymbol{\varepsilon}_t
\tag{5}
$$

단, $$\Delta \boldsymbol{X}_t$$는 변수들의 1차 차분 항의 ($$k \times 1$$)벡터, $$\boldsymbol{\Gamma}_i$$와 $$\boldsymbol{\Pi}$$는 ($$k \times k$$) 계수행렬

#### 장기균형 행렬

&#8203;$$rank(\boldsymbol{\Pi})$$는 (5)를 이해하고 공적분검증에 있어 필수적이다. $$\boldsymbol{\Pi}$$가 (4)의 VAR모형 변수들의 장기균형관계에 대한 정보를 포함하고 있기 때문에, $$\boldsymbol{\Pi}$$는 흔히 장기균형행렬이라고도 한다.

#### 공적분 벡터의 존재성과 장기균형행렬의 랭크

다음의 VECM 모형을 살펴보자.

$$
\Delta \boldsymbol{X}_t = \boldsymbol{\Pi} \boldsymbol{X}_{t-1} + \sum_{i=1}^{p-1} \boldsymbol{\Gamma}_i \Delta \boldsymbol{X}_{t-i} + \boldsymbol{\varepsilon}_t
$$

&#8203;$$\Delta \boldsymbol{X}_t$$와 오차항 $$\boldsymbol{\varepsilon}_t$$는 모두 정상과정($$I(0)$$)이다. 따라서 모형의 양변이 균형을 이루기 위해서는 우변의 $$\boldsymbol{\Pi}\boldsymbol{X}_{t-1}$$ 역시 정상과정이어야만 한다. 변수 벡터 $$\boldsymbol{X}_{t-1}$$이 불안정적인 1차 적분 변수($$I(1)$$)들로 구성되어 있으므로, $$\boldsymbol{\Pi}\boldsymbol{X}_{t-1}$$이 정상적이 되기 위한 조건은 장기균형행렬 $$\boldsymbol{\Pi}$$가 선형적으로 독립적인 행을 얼마나 가지느냐, 즉 행렬의 rank에 의해 결정된다.

내생변수의 수를 $$k$$라 하고, 행렬 $$\boldsymbol{\Pi}$$의 랭크를 $$r$$이라 할 때 ($$r = \text{rank}(\boldsymbol{\Pi})$$), 공적분 벡터의 존재성은 다음 세 가지 경우로 분류할 수 있다.

**1. $$r = 0$$ 인 경우**

&#8203;$$\boldsymbol{\Pi}$$가 영행렬($$\boldsymbol{\Pi} = \boldsymbol{0}$$)임을 의미한다. 변수들 간에 선형결합을 통해 정상성을 확보할 수 있는 공적분 관계가 전혀 존재하지 않는다. 이 경우 공적분 벡터는 존재하지 않으며, 모형은 단순히 1차 차분된 변수들의 다변량 모형인 VAR($$p-1$$) 모형으로 축소된다.

**2. $$0 < r < k$$ 인 경우**

일반적인 공적분 모형의 형태다. 불안정적인 $$k$$개의 변수들 사이에 $$r$$개의 선형독립인 공적분 벡터가 존재한다. 이때 행렬 $$\boldsymbol{\Pi}$$는 $$\boldsymbol{\Pi} = \boldsymbol{\alpha}\boldsymbol{\beta}'$$의 형태로 분해(Granger 표현의 정리)될 수 있으며, 행렬 $$\boldsymbol{\beta}$$를 구성하는 열벡터들이 바로 $$r$$개의 공적분 벡터가 된다.

**3. $$r = k$$ 인 경우 (full rank)**

&#8203;$$\boldsymbol{\Pi}$$가 non-singular matrix일 때다. 이 경우 $$\boldsymbol{X}_{t-1}$$을 구성하는 모든 변수 자체가 이미 정상과정($$I(0)$$)이라는 것을 내포한다. 변수들이 원래부터 안정적이므로 공적분 분석은 필요하지 않으며, 차분할 필요 없이 수준변수(level variables) 그대로를 이용해 일반적인 VAR 모형을 추정해야 한다.

#### Granger 표현 정리

Granger 표현 정리에 의하면, $$rank(\boldsymbol{\Pi})$$가 적분된 변수들의 수보다 작다면(즉, $$r \leq k-1$$), $$\boldsymbol{\Pi} = \boldsymbol{\alpha}\boldsymbol{\beta}^T$$이고 $$\boldsymbol{\beta}^T\boldsymbol{X}_t$$의 모든 성분들이 정상적이 되도록 하는 ($$k \times r$$)행렬인 $$\boldsymbol{\alpha}$$와 $$\boldsymbol{\beta}$$가 존재한다.

$$
\Delta \boldsymbol{X}_t = \underbrace{\boldsymbol{\alpha}\boldsymbol{\beta}^T\boldsymbol{X}_{t-1}}_{\boldsymbol{\alpha} \boldsymbol{e}_{t-1}} + \sum_{i=1}^{p-1} \boldsymbol{\Gamma}_i \Delta \boldsymbol{X}_{t-i} + \boldsymbol{\varepsilon}_t
\tag{6}
$$

행렬 $$\boldsymbol{\beta}$$는 공적분벡터들로 이루어진 행렬이다. 행렬 $$\boldsymbol{\alpha}$$의 성분들은 공적분벡터의 성분에 대한 가중치로 해석될 수 있다.

### VECM의 구성요소

(6)은 장기균형관계를 회복하기 위해 적분된 변수들이 보여주는 단기적·동태적 특성을 나타내기에 **벡터오차수정모형**이라고 한다.

- 공적분관계는 $$\boldsymbol{X}_t$$를 이루는 적분된 변수들이 장기균형관계를 갖는다는 것을 의미한다. $$\boldsymbol{X}_t$$의 모든 변수들이 언제나 균형상태를 유지하고 있다면, $$\boldsymbol{e}_{t-1}$$값은 모두 0이다.
- 만약 무작위 충격으로 장기균형으로부터의 일시적인 이탈이 발생하면, 즉 $$\boldsymbol{e}_{t-1}$$이 0이 아니라면, 오차수정 메커니즘이 균형을 회복하기 위해 작동한다.
- $$\boldsymbol{\alpha}$$의 계수들은 그러한 오차수정의 속도를 나타내주는 것으로서, 수정속도계수 또는 수정계수라고 부른다. $$\boldsymbol{\alpha}$$의 계수들 중 큰 값은 그에 대응하는 변수가 균형값을 회복하기 위해 더 빨리 변화한다는 것을 의미한다.
- 만약 $$\boldsymbol{\alpha}$$의 계수 중 어떤 것이 0에 가깝다면, 이에 대응하는 변수는 $$t-1$$시점에 발생한 장기균형으로부터의 이탈에 대해 거의 반응하지 않는다는 것을 의미한다.

## VECM 해석

다음과 같은 벡터오차수정모형이 주어져 있다고 가정하자.

$$
\begin{pmatrix}
\Delta x_t \\
\Delta y_t
\end{pmatrix}
=
\begin{pmatrix}
\alpha_1 \\
\alpha_2
\end{pmatrix}
(x_{t-1} + \beta y_{t-1} + c)
+
\begin{pmatrix}
\phi_{11,1} & \phi_{12,1} \\
\phi_{21,1} & \phi_{22,1}
\end{pmatrix}
\begin{pmatrix}
\Delta x_{t-1} \\
\Delta y_{t-1}
\end{pmatrix}
+
\begin{pmatrix}
\phi_{11,2} & \phi_{12,2} \\
\phi_{21,2} & \phi_{22,2}
\end{pmatrix}
\begin{pmatrix}
\Delta x_{t-2} \\
\Delta y_{t-2}
\end{pmatrix}
+
\begin{pmatrix}
\varepsilon_{x,t} \\
\varepsilon_{y,t}
\end{pmatrix}
\tag{7}
$$

#### 선형결합의 정상성

정규화된 공적분 방정식은 $$e_t = x_{t-1} + \beta y_{t-1} + c$$이다. 즉, 공적분벡터 $$(1, \beta)$$를 갖는 $$x_t$$와 $$y_t$$의 선형결합은 평균이 $$c$$인 정상과정을 갖는다.

#### 조정계수

조정계수는 공적분된 변수들이 장기균형을 회복할 때의 조정속도를 보여준다. 조정계수의 의미를 좀 더 명확하게 하기 위해서, (7)을 다음과 같이 다시 쓰자.

$$
\begin{pmatrix}
\Delta x_t \\
\Delta y_t
\end{pmatrix}
=
\begin{pmatrix}
\alpha_1 e_{t-1} \\
\alpha_2 e_{t-1}
\end{pmatrix}
+
\begin{pmatrix}
\phi_{11,1} & \phi_{12,1} \\
\phi_{21,1} & \phi_{22,1}
\end{pmatrix}
\begin{pmatrix}
\Delta x_{t-1} \\
\Delta y_{t-1}
\end{pmatrix}
+
\begin{pmatrix}
\phi_{11,2} & \phi_{12,2} \\
\phi_{21,2} & \phi_{22,2}
\end{pmatrix}
\begin{pmatrix}
\Delta x_{t-2} \\
\Delta y_{t-2}
\end{pmatrix}
+
\begin{pmatrix}
\varepsilon_{x,t} \\
\varepsilon_{y,t}
\end{pmatrix}
\tag{8}
$$

시점 $$t-1$$에 $$x_{t-1}$$과 $$y_{t-1}$$의 값이 장기균형으로부터 이탈한다면, 첫째 항에 있는 $$e_{t-1}$$은 0이 아닌 값을 갖게 된다. 균형으로부터 이탈이 있을 때, 한 단위의 $$e_{t-1}$$은 평균적으로 $$x_t$$가 $$\alpha_1$$만큼 변하고, $$y_t$$가 $$\alpha_2$$만큼 변하도록 만든다.

#### 그랜저 인과성

(7)의 둘째 및 셋째 항에 있는 $$\phi$$값들은 $$x$$와 $$y$$ 사이의 그랜저 인과성과 선도/지연 관계에 대해 몇 가지 시사점을 갖고 있다. (7)과 (8)을 정리하면 다음 식이 얻어진다.

$$
\begin{aligned}
\Delta x_t &= \alpha_1 e_{t-1} + \phi_{11,1} \Delta x_{t-1} + \phi_{11,2} \Delta x_{t-2} + \phi_{12,1} \Delta y_{t-1} + \phi_{12,2} \Delta y_{t-2} + \varepsilon_{x,t} \qquad \text{(9)} \\
\Delta y_t &= \alpha_2 e_{t-1} + \phi_{21,1} \Delta x_{t-1} + \phi_{21,2} \Delta x_{t-2} + \phi_{22,1} \Delta y_{t-1} + \phi_{22,2} \Delta y_{t-2} + \varepsilon_{y,t} \qquad \text{(10)}
\end{aligned}
$$

- (9)에서 $$\phi_{12,1}$$과 $$\phi_{12,2}$$가 0과 다르지 않으면, $$y$$는 $$x$$를 그랜저 인과하지 않는다.
- (10)에서 $$\phi_{21,1}$$과 $$\phi_{21,2}$$가 0과 다르지 않으면, $$x$$는 $$y$$를 그랜저 인과하지 않는다.
- (9)에서 $$\phi_{12,1}$$이 유의한 값이고 $$\phi_{12,2}$$가 유의하지 않다면, $$y$$가 $$x$$를 1기간만큼 선도한다고 해석될 수 있다.
- $$\phi_{12,2}$$가 유의하다면, $$y$$가 2기간만큼 $$x$$를 선도한다는 것을 의미한다.
