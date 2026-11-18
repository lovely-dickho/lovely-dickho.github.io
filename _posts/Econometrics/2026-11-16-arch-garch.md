---
title: "[Econometrics] ARCH and GARCH"
date: 2026-11-16 18:55:00 +0900
categories: [Econometrics, Notes]
tags: [ARCH, GARCH, 변동성, TARCH, EGARCH]
math: true
toc: true
---

## 총분산의 법칙 (Law of Total Variance)

확률변수 $$Y$$의 무조건부 분산은 그룹 내 분산의 기댓값과 그룹 간 분산의 합으로 분해된다.

$$
Var(Y) = E[Var(Y|X)] + Var(E[Y|X])
$$

시계열 오차항 $$u_t$$의 경우, 조건부 기댓값이 $$E[u_t | \Omega_{t-1}] = 0$$이므로 두 번째 항이 0이 된다. 따라서 무조건부 분산은 조건부 분산의 기댓값으로 단순화된다.

$$
Var(u_t) = E[Var(u_t | \Omega_{t-1})] + Var(0) = E[h_t]
$$

## ARCH

ARCH모형은 $$t$$시점의 조건부분산이 잔차제곱의 시차값에 의존하는 조건부변동성모형이다. ARCH(1)은 다음 식과 같이 $$t$$기의 조건부분산이 $$t-1$$기의 잔차제곱에 의존하는 모형이다.

$$
\begin{aligned}
\text{평균방정식} &: Y_t = \mu_t + u_t, \text{ 단 } u_t = \sqrt{h_t} \varepsilon_t \qquad \text{(1)} \\
\text{분산방정식} &: h_t = \alpha_0 + \alpha_1 u_{t-1}^2 \qquad \text{(2)}
\end{aligned}
$$

단, $$\mu_t$$는 조건부평균. $$u_t$$는 잔차항. $$\varepsilon_t$$는 평균이 0이고 분산이 1인 백색잡음. $$h_t$$는 $$u_t$$의 조건부분산

- 평균방정식은 $$Y_t$$의 확률과정을 정의한다. 조건부평균 $$\mu_t$$는 $$Y_t$$의 기대치를 나타낸다. 평균방정식은 몇 개의 독립변수를 포함하고 있는 회귀모형이거나 정상적 ARMA모형일 수 있다.

  - $$Y_t$$가 ARMA(1,1)과정을 갖는다면, 평균방정식의 구체적인 형태는 $$Y_t = \phi Y_{t-1} + u_t + \theta u_{t-1}$$이 될 것이고, 이런 모형을 ARMA(1,1)-ARCH(1)모형이라고 한다
- ARCH에서는 $$u_t$$가 시간가변적 분산 $$h_t$$의 제곱근과 평균이 0이고 분산이 1인 백색잡음 $$\varepsilon_t$$의 곱이기 때문에, $$u_t$$는 0의 평균과 시간가변적 분산 $$h_t$$를 갖는다.
- 분산방정식 (2)는 시간가변적 분산 $$h_t$$를 정의하는 식이다. ARCH(1) 모형에서 $$h_t$$는 잔차제곱의 시차값에 의해 결정된다.
- 분산은 음(-)일 수 없기 때문에, 모수 $$\alpha_0$$과 $$\alpha_1$$ 역시 음(-)일 수 없다. 분산방정식의 식 (2)는 $$t-1$$기에 커다란 충격, 즉 큰 $$u_{t-1}^2$$이 주어지면, $$t$$기에도 커다란 변동성이 예상된다는 것을 시사하고 있다.

## ARCH(1)모형

오차항 $$u_t$$가 다음과 같은 구조를 가진다고 가정한다.

$$
u_t = \sqrt{h_t} \varepsilon_t
$$

여기서 $$\varepsilon_t \sim i.i.d. N(0, 1)$$이며, 과거 정보집합을 $$\Omega_{t-1}$$이라 할 때 조건부 분산 $$h_t$$은 다음과 같다.

$$
h_t = Var(u_t | \Omega_{t-1}) = \alpha_0 + \alpha_1 u_{t-1}^2
$$

조건부 분산이 양수이고 정상성을 만족하기 위해 $$\alpha_0 > 0$$, $$0 \le \alpha_1 < 1$$을 가정한다.

#### 무조건부 평균 도출

반복기댓값의 법칙을 이용하여 도출한다.

$$
E[u_t] = E[E[u_t | \Omega_{t-1}]] = E[\sqrt{h_t} E[\varepsilon_t | \Omega_{t-1}]] = E[\sqrt{h_t} \cdot 0] = 0
$$

#### 무조건부 분산

무조건부 기댓값이 0이므로 분산은 $$E[u_t^2]$$와 같으며, 앞서 살펴본 총분산의 법칙에 의해서도 $$Var(u_t) = E[h_t]$$로 단순화된다.

$$
\begin{aligned}
E[u_t^2] &= E[E[u_t^2 | \Omega_{t-1}]] = E[h_t E[\varepsilon_t^2 | \Omega_{t-1}]] = E[h_t] \\
E[h_t] &= E[\alpha_0 + \alpha_1 u_{t-1}^2] = \alpha_0 + \alpha_1 E[u_{t-1}^2]
\end{aligned}
$$

공분산 정상성에 의해 $$E[u_t^2] = E[u_{t-1}^2] = \sigma^2$$이 성립하므로,

$$
\sigma^2 = \alpha_0 + \alpha_1 \sigma^2 \implies Var(u_t) = \frac{\alpha_0}{1 - \alpha_1}
$$

#### 무조건부 자기공분산

시차 $$k > 0$$에 대하여,

$$
Cov(u_t, u_{t-k}) = E[u_t u_{t-k}] = E[E[u_t u_{t-k} | \Omega_{t-1}]] = E[u_{t-k} E[u_t | \Omega_{t-1}]] = 0
$$

따라서 오차항은 무상관 백색잡음 형태를 띤다.

#### 3차 적률과 대칭성 (Skewness)

분포의 비대칭성을 확인하기 위해 3차 적률 $$E[u_t^3]$$을 도출한다. 마찬가지로 반복기댓값의 법칙을 적용한다.

$$
E[u_t^3] = E[E[u_t^3 | \Omega_{t-1}]] = E[h_t^{3/2} E[\varepsilon_t^3 | \Omega_{t-1}]]
$$

&#8203;$$\varepsilon_t$$는 표준정규분포를 따르므로 평균을 중심으로 완벽히 대칭이며, 3차 적률 $$E[\varepsilon_t^3] = 0$$이다.

$$
E[u_t^3] = E[h_t^{3/2} \cdot 0] = 0
$$

왜도는 $$E[u_t^3] / (E[u_t^2])^{3/2} = 0$$이 되어, ARCH 모형의 오차항은 무조건부로도 대칭분포를 가짐을 알 수 있다.

#### 4차 적률과 두터운 꼬리(Fat tail) 증명

금융 시계열은 정규분포보다 꼬리가 두꺼운 현상(첨도 $$> 3$$)을 자주 보인다. ARCH(1) 모형의 첨도(Kurtosis)를 도출하여 이를 보이자.

첨도 $$K$$는 다음과 같이 정의된다.

$$
K = \frac{E[u_t^4]}{(E[u_t^2])^2}
$$

우선 4차 적률 $$m_4 = E[u_t^4]$$를 구한다. $$\varepsilon_t \sim N(0, 1)$$이므로 $$E[\varepsilon_t^4] = 3$$이다.

$$
\begin{aligned}
E[u_t^4 | \Omega_{t-1}] &= E[\varepsilon_t^4 h_t^2 | \Omega_{t-1}] = 3h_t^2 \\
m_4 &= E[E[u_t^4 | \Omega_{t-1}]] = 3E[h_t^2]
\end{aligned}
$$

여기서 $$E[h_t^2]$$를 전개한다.

$$
E[h_t^2] = E[(\alpha_0 + \alpha_1 u_{t-1}^2)^2] = \alpha_0^2 + 2\alpha_0\alpha_1 E[u_{t-1}^2] + \alpha_1^2 E[u_{t-1}^4]
$$

정상성에 의해 $$E[u_{t-1}^2] = \frac{\alpha_0}{1-\alpha_1}$$이고 $$E[u_{t-1}^4] = m_4$$이므로 대입한다.

$$
E[h_t^2] = \alpha_0^2 + 2\alpha_0\alpha_1 \left( \frac{\alpha_0}{1 - \alpha_1} \right) + \alpha_1^2 m_4
$$

이를 $$m_4 = 3E[h_t^2]$$ 식에 대입하여 $$m_4$$에 대해 정리한다.

$$
\begin{aligned}
m_4 &= 3 \left( \alpha_0^2 + \frac{2\alpha_0^2\alpha_1}{1 - \alpha_1} + \alpha_1^2 m_4 \right) \\
m_4 (1 - 3\alpha_1^2) &= 3\alpha_0^2 \left( 1 + \frac{2\alpha_1}{1 - \alpha_1} \right) = 3\alpha_0^2 \left( \frac{1 + \alpha_1}{1 - \alpha_1} \right) \\
m_4 &= \frac{3\alpha_0^2 (1 + \alpha_1)}{(1 - \alpha_1)(1 - 3\alpha_1^2)}
\end{aligned}
$$

(단, 4차 적률이 존재하기 위해 $$3\alpha_1^2 < 1$$ 조건이 추가로 요구된다.)

이제 첨도 $$K$$를 계산한다. $$(E[u_t^2])^2 = \left( \frac{\alpha_0}{1 - \alpha_1} \right)^2 = \frac{\alpha_0^2}{(1 - \alpha_1)^2}$$ 이다.

$$
K = \frac{ \frac{3\alpha_0^2 (1 + \alpha_1)}{(1 - \alpha_1)(1 - 3\alpha_1^2)} }{ \frac{\alpha_0^2}{(1 - \alpha_1)^2} } = \frac{3(1 + \alpha_1)(1 - \alpha_1)}{1 - 3\alpha_1^2} = \frac{3(1 - \alpha_1^2)}{1 - 3\alpha_1^2}
$$

&#8203;$$\alpha_1 > 0$$ 이므로, 분자의 $$1 - \alpha_1^2$$는 분모의 $$1 - 3\alpha_1^2$$보다 항상 크다. 따라서 다음 부등식이 성립한다.

$$
K = 3 \left( \frac{1 - \alpha_1^2}{1 - 3\alpha_1^2} \right) > 3
$$

결과적으로 ARCH 모형은 오차항 $$u_t$$의 조건부분포가 정규분포를 따르더라도 3차 적률이 0이므로 대칭분포를 가지지만, 변동성 군집 현상으로 인해 무조건부분포는 정규분포보다 꼬리가 두터운(Fat tail) 특성을 지님을 알 수 있다.

### ARCH(1)과 AR(1)의 관계

ARCH(1) 모형에서 오차항의 제곱인 $$u_t^2$$은 1차 자기회귀모형인 AR(1) 과정을 따른다. $$u_t$$의 조건부분산 $$h_t$$는 다음과 같이 정의된다.

$$
h_t = E_{t-1}(u_t^2) = \alpha_0 + \alpha_1 u_{t-1}^2
$$

여기서 새로운 오차항 $$v_t$$를 실제 오차항의 제곱($$u_t^2$$)과 조건부분산($$h_t$$)의 차이로 정의하자.

$$
v_t = u_t^2 - h_t
$$

이때 $$v_t$$의 조건부 기댓값은 $$E_{t-1}(v_t) = E_{t-1}(u_t^2) - h_t = h_t - h_t = 0$$이 되어 백색잡음(White Noise)의 성질을 만족한다. 위 식을 $$u_t^2$$에 대해 정리하고 $$h_t$$를 대입하면 다음과 같다.

$$
u_t^2 = h_t + v_t = \alpha_0 + \alpha_1 u_{t-1}^2 + v_t
$$

결과적으로 $$u_t^2$$은 자기 자신의 과거 값인 $$u_{t-1}^2$$에 의존하는 전형적인 AR(1) 모형의 형태를 띠게 된다.

### ARCH(q)모형

ARCH($$q$$)모형에서 조건부분산은 잔차제곱의 $$q$$번째까지의 시차값에 의해 결정된다.

$$
h_t = \alpha_0 + \alpha_1 u_{t-1}^2 + \cdots + \alpha_q u_{t-q}^2 \quad where \quad \alpha_0, \alpha_1, \cdots, \alpha_q >0
\tag{3}
$$

어떤 경우에는 차수 $$q$$가 너무 커서 간결성의 원칙을 지키지 못할 수도 있다.

&#8203;$$Y_t$$의 변동성에 중요한 영향을 미치는 외생변수들이 있다면, 그런 변수도 분산방정식에 포함될 수 있다. 만일 $$X_t$$가 그런 외생변수라면, ARCH($$q$$)모형은 다음과 같이 표현된다.

$$
h_t = \alpha_0 + \lambda X_t + \alpha_1 u_{t-1}^2 + \cdots + \alpha_q u_{t-q}^2
$$

## ARCH Topics

평균방정식의 구체적인 형태가 정해져 있고 잔차항이 주어져 있으면, **ARCH($$q$$)모형의 차수 결정, 모형의 모수 추정, 모형설정의 검증**의 3단계를 거쳐 ARCH모형을 추정할 수 있다.

### 모형선택

ARCH모형을 추정하는 첫 번째 단계는 평균방정식의 잔차항 $$u_t$$에 ARCH효과가 존재하는지와 ARCH($$q$$)모형의 차수가 무엇인지를 확인하는 일이다.

#### ARCH 효과 검증

- **$$u_t^2$$ 시계열의 Ljung-Box $$Q$$통계량**: 대부분의 경우 $$Q$$통계량은 ARCH모형의 차수에 대해서 분명한 신호를 주지 못함
- **LM test**: 식 (4)의 회귀식에서 $$H_0 : \alpha_1 = \cdots = \alpha_q = 0$$을 검증

  $$
  u_t^2 = \alpha_0 + \alpha_1 u_{t-1}^2 + \cdots + \alpha_q u_{t-q}^2 + \eta_t
  \tag{4}
  $$

  단, $$\eta_t$$는 오차항임.

귀무가설의 검증을 위해 두 가지 통계량이 사용된다: $$F$$통계량, 자유도가 $$q$$인 카이제곱분포를 갖는 $$TR^2$$ 통계량

**①** 어떤 통계량을 사용하든 일단 충분히 큰 $$q$$값을 선택하여 회귀분석하여야 한다.

**②** 만약 $$F$$나 $$TR^2$$의 값이 귀무가설을 기각하고 모든 $$\alpha$$값이 0보다 유의하게 크다면, $$u_t$$는 차수가 $$q$$인 ARCH효과를 갖고 있다고 결론내린다.

**③** $$F$$나 $$TR^2$$ 값이 귀무가설을 기각하지만 마지막의 몇 개의 $$\alpha$$값들만 0보다 유의하게 크지 않다면, 더 적은 차수의 회귀식을 추정해야 한다.

### 추정

ARCH모형은 MLE로 추정한다. 만일 $$u_t$$가 평균이 0이고 분산이 $$h_t$$인 정규분포를 갖는다고 가정하면, 우도함수는 다음과 같다.

$$
L(\alpha_0, \cdots, \alpha_q) = \prod_{t=q+1}^{T} \frac{1}{\sqrt{2\pi h_t}} e^{-0.5[u_t^2 / h_t]} \times f(u_1, \cdots, u_q | \alpha_0, \cdots, \alpha_q)
\tag{5}
$$

여기서 $$f(u_1, \cdots, u_q | \alpha_0, \cdots, \alpha_q)$$는 모수 $$\alpha_0, \cdots, \alpha_q$$의 값이 주어져 있을 때의 $$\{u_1, \cdots, u_q\}$$의 결합확률이다.

시점 $$t$$마다 조건부분산을 추정해야 하므로, 위의 우도함수에서는 일정한 분산 $$\sigma^2$$ 대신 $$h_t$$가 사용된다. 실제 자료를 사용할 경우, $$h_t$$는 실제의 분산방정식으로, $$u_t$$는 평균방정식의 잔차로 대체된다.

로그우도함수는 다음과 같다.

$$
\log L = \sum_{t=q+1}^{T} \left[ -\frac{1}{2} \ln(2\pi) - \frac{1}{2} \ln(h_t) - \frac{1}{2} \frac{u_t^2}{h_t} \right] + \ln [f(u_1, \cdots, u_q)]
\tag{6}
$$

실제의 추정에서 $$\ln[f(u_1, \cdots, u_q)]$$의 값을 계산하기 어렵기 때문에, 관찰치 수 $$T$$가 충분히 클 때에는 이 부분을 무시하는 경우가 많다.

### 모형의 진단

- 식 (1)에서 잔차항 $$u_t$$는 $$\sqrt{h_t} \varepsilon_t$$로 정의했다.
- 이것은 ARCH모형이 $$h_t$$를 적절하게 추정해 준다면, 표준화된 시계열 $$\hat{\varepsilon}_t = u_t / \sqrt{h_t}$$는 백색잡음이어야 한다.
- $$\hat{\varepsilon}_t$$의 시계열이 백색잡음인지를 검증함으로써, 추정된 ARCH모형의 적절성을 확인할 수 있다.

## GARCH(1,1) 모형으로의 일반화

ARCH 모형은 긴 시차의 변동성을 포착하기 위해 많은 파라미터가 필요하다는 단점이 있다. 이를 해결하기 위해 과거의 조건부 분산 자체를 설명변수로 포함시킨 모형이 GARCH 모형이다. GARCH(1,1)의 조건부 분산 식은 다음과 같다.

$$
h_t = \omega + \alpha u_{t-1}^2 + \beta h_{t-1}
$$

무조건부 분산을 구하기 위해 양변에 기댓값을 취한다.

$$
E[h_t] = \omega + \alpha E[u_{t-1}^2] + \beta E[h_{t-1}]
$$

정상성 가정 하에 $$E[u_t^2] = E[h_t] = \sigma^2$$이 성립하므로,

$$
\sigma^2 = \omega + \alpha \sigma^2 + \beta \sigma^2 \implies Var(u_t) = \frac{\omega}{1 - \alpha - \beta}
$$

정상성을 만족하기 위한 조건은 $$\alpha + \beta < 1$$이다.

### 모형

가장 단순하고 가장 널리 이용되는 GARCH모형은 식 (7)와 식 (8)와 같은 형태의 GARCH(1,1)모형이다.

$$
\begin{aligned}
\text{평균방정식} &: Y_t = \mu_t + u_t, \quad \text{단, } u_t = \sqrt{h_t} \varepsilon_t \qquad \text{(7)} \\
\text{분산방정식} &: h_t = \omega + \alpha u_{t-1}^2 + \beta h_{t-1} \qquad \text{(8)}
\end{aligned}
$$

#### 비음제약

ARCH모형과 유사한 제약조건들이 GARCH(1,1)모형에도 적용된다.

- $$h_t$$가 비음이어야 하므로, 모수 $$\omega, \alpha$$ 및 $$\beta$$가 음(-)이 아니어야 한다.
- $$u_t$$가 일정한 무조건부 분산을 갖기 위해서, $$\beta \leq 1$$와 $$(\alpha + \beta) < 1$$의 조건을 만족시켜야 한다.

#### 일반화

GARCH(1,1)모형은 다음과 같은 GARCH(q,p)모형으로 일반화될 수 있다.

$$
\text{분산방정식} : h_t = \omega + \sum_{i=1}^{q} \alpha_i u_{t-i}^2 + \sum_{i=1}^{p} \beta_i h_{t-i}
\tag{9}
$$

이 모형은 잔차제곱 항인 $$u_t^2$$의 ARMA모형과 동일한 형식이다. GARCH(1,1)모형에서와 마찬가지로 다음 조건을 만족해야 한다.

- $$h_t$$가 비음이기 위해서 모수 $$\omega, \alpha_i, \cdots, \beta_i$$ 등이 비음이어야한다.
- $$u_t$$가 일정한 무조건부분산을 갖기 위해서 $$\sum_{i=1}^{q} \alpha_i + \sum_{i=1}^{p} \beta_i < 1$$이어야 한다는 제약조건이 충족되어야 한다.

### I-GARCH

ARCH($$q$$)와 GARCH($$q,p$$)모형의 분산방정식은 $$u_t^2$$에 관한 ARMA모형과 비슷하다. 이렇게 표현된 식의 AR항들이 단위근을 가질 경우, 무작위적 충격이 주어지면 그것의 제곱값이 시계열과정에 영속적인 영향을 미치게 된다.

GARCH(1,1)모형의 분산방정식에서 $$a_1 + b_1$$이 1과 같으면, 이 모형은 적분된 GARCH(1,1)이라고 하며, 이 모형의 분산방정식은 다음과 같이 나타내진다.

$$
h_t = \omega + (1 - \beta) u_{t-1}^2 + \beta h_{t-1}
\tag{10}
$$

적분된 GARCH모형에서는 $$u_t$$ 및 $$Y_t$$의 무조건부분산이 유한한 값을 갖지 않는다는 사실은 쉽게 증명될 수 있다.

**증명.** $$\eta_t$$를 $$u_{t}^2 - h_t$$로 정의하자. GARCH모형은 다음 식으로 다시 쓸 수 있다.

$$
u_t^2 = \omega + \sum_{i=1}^{\max(q,p)} (\alpha_i + \beta_i) u_{t-i}^2 + \eta_t - \sum_{i=1}^{p} \beta_i \eta_{t-i}
$$

이 식은 $$\max(q,p)$$개의 AR 항과 $$p$$개의 MA 항으로 이루어진 $$u_t^2$$의 ARMA모형이다. 여기서 조건부분산은 $$E(u_t^2) = \omega / [1 - \sum_{i=1}^{\max(q,p)} (\alpha_i + \beta_i)]$$으로 도출된다(단, 분모는 양의 값을 가져야 한다). 따라서 무조건부분산이 음이 아니고 유한한 값을 갖기 위한 제약조건을 이 무조건부분산의 식에서 도출할 수 있다.

위 식으로부터 $$E(u_t^2) = \omega / [1 - (\alpha + \beta)]$$이 얻어진다. 이 식에서 $$\alpha + \beta$$이 1이면, 무조건부분산 $$E(u_t^2)$$은 무한한 값을 갖는다. ∎

그 결과로서 $$u_t^2$$의 예상치 못한 충격이 발생하면, 이 충격은 시계열과정에 대해 영속적인 영향을 미치기 때문에, 조건부분산이 새로운 수준으로 이동하도록 만들 수 있다.

## ARCH-M and GARCH-M

재무학에서 가장 중요한 원칙 중의 하나는 위험·수익의 교환관계이다. 사람들은 더 많은 위험을 부담함으로써 더 높은 수익률이 얻어지기를 기대한다. 따라서 시간가변적 위험이 시간가변적 기대수익률을 가져올 것이라는 것은 충분히 합리적인 생각이라고 할 수 있다.

이 아이디어는 Engle, Lilien & Robins(1987)의 ARCH-M모형과 French, Schwert & Stambaugh(1987)의 GARCH-M모형으로 구체화되었다. 이들의 모형 이름에서 M은 평균을 의미하는 것으로서, 평균방정식이 조건부변동성을 독립변수로 포함하고 있음을 나타낸다.

예를 들어, GARCH(1,1)-M모형은 다음과 같이 표현된다.

$$
\begin{aligned}
\text{평균방정식} &: Y_t = \mu_t + \gamma h_t + u_t \qquad \text{(11)} \\
\text{분산방정식} &: h_t = \omega + \alpha u_{t-1}^2 + \beta h_{t-1} \qquad \text{(12)}
\end{aligned}
$$

위 식에서 조건부분산 $$h_t$$가 평균방정식에 포함되어 있고 계수 $$\gamma$$는 양의 값을 가질 것으로 기대된다. 조건부표준편차 $$\sqrt{h_t}$$가 조건부분산 대신 사용될 수도 있다.

## TARCH

실증적 결과는 주가가 하락하면 뒤이어 더 큰 변동성이 나타난다는 것을 보여준다. 변동성 ‘레버리지 효과’에 의하면, 주식의 시장가격이 하락하면, 재무레버리지 - 즉, 부채와 자기자본의 시장가치의 비율로 계산되는 부채비율-가 증가하고, 이에 따라 자기자본 소유자들의 재무위험이 증가하기 때문에 이런 현상이 나타난다고 한다.

이런 현상을 모형화하기 위해 Glosten, Jaganathan & Runkle(1993)은 다음과 같은 TARCH모형을 제시하였다.

$$
h_t = \omega + \alpha u_{t-1}^2 + \gamma u_{t-1}^2 d_{t-1} + \beta h_{t-1}
\tag{13}
$$

단, $$d_{t-1}$$은 나쁜 뉴스($$u_{t-1} < 0$$)일 때 $$d_{t-1} = 1$$이고, 좋은 뉴스($$u_{t-1} \geq 0$$)일 때 $$d_{t-1} = 0$$인 더미변수

더미변수 값을 정의하는 기준점은 식 (13)에서 0이다. 나쁜 뉴스가 발생하면 1단위의 $$u_{t-1}^2$$에 대한 변동성의 반응은 $$\alpha + \gamma$$이고, 좋은 뉴스가 일어나면 $$\alpha$$이다. 따라서 TARCH모형은 변동성의 비대칭적 반응을 고려하고 있지만, 음의 추정치가 얻어질 가능성은 배제하지 않고 있다.

## EGARCH

Nelson(1991)은 EGARCH모형이라고 부르는 다른 접근방법을 제안하였다. EGARCH(1,1)모형의 형태는 다음과 같다.

$$
\begin{aligned}
\text{평균방정식} &: Y_t = \mu_t + u_t, \quad \text{단, } u_t = \sqrt{h_t} \varepsilon_t \qquad \text{(14)} \\
\text{분산방정식} &: \ln(h_t) = \omega + \alpha g(\varepsilon_{t-1}) + \beta \ln(h_{t-1}) \qquad \text{(15)}
\end{aligned}
$$

단, $$g(\varepsilon_{t-1}) = \theta \varepsilon_{t-1} + \gamma [|\varepsilon_{t-1}| - E(|\varepsilon_{t-1}|)]$$

분산방정식의 좌변이 $$h_t$$의 로그값이므로, 로그를 풀면 $$h_t = e^{(\text{우변})}$$이 된다.

- 우변의 값이 음인 경우라 하더라도 $$e^{(\text{우변})}$$의 값은 언제나 양수가 되기 때문에 조건부분산 $$h_t$$ 역시 언제나 양의 값을 갖는다.
- 우변이 음의 값을 가져도 괜찮기 때문에 분산방정식에서 모수들이 음의 값을 가져도 문제될 것이 없다.

&#8203;$$g(\varepsilon_{t-1})$$은 GARCH모형에서 $$u_{t-1}^2$$의 역할을 한다.

- $$\varepsilon_{t-1}$$은 i.i.d.이며 0의 평균과 1의 분산을 갖는 백색잡음으로서, $$\varepsilon_{t-1} = u_{t-1}/\sqrt{h_{t-1}}$$에 의해 얻어진 백색잡음이다.
- $$\varepsilon_{t-1}$$의 절댓값인 $$|\varepsilon_{t-1}|$$ 역시 i.i.d.인 백색잡음과정이다.
- $$\varepsilon_{t-1}$$이 정규분포를 갖는다면, 절댓값의 기대치는 $$E(|\varepsilon_{t-1}|) = \sqrt{2/\pi}$$이다.

절댓값 연산기호는 함수 $$g$$가 조건부분산의 비대칭적 반응을 나타낼 수 있도록 만든다. $$g$$를 $$\varepsilon_{t-1}$$이 0보다 큰 경우와 작은 경우로 나누어 전개해 보자.

$$
g(\varepsilon_{t-1}) =
\begin{cases}
(\theta + \gamma)\varepsilon_{t-1} - \gamma E(|\varepsilon_{t-1}|) & \varepsilon_{t-1} \geq 0 \text{이면} \\
(\theta - \gamma)\varepsilon_{t-1} - \gamma E(|\varepsilon_{t-1}|) & \varepsilon_{t-1} < 0 \text{이면}
\end{cases}
\tag{16}
$$

&#8203;$$\theta$$가 음수이면 전날의 충격이 음수일 때 $$g(\varepsilon_{t-1})$$의 값과 조건부분산이 증가하는 경향이 있다.
