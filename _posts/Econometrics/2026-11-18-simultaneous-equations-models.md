---
title: "[Econometrics] Simultaneous Equations Models"
date: 2026-11-18 12:25:00 +0900
categories: [Econometrics, Notes]
tags: [연립방정식, 식별, 2SLS, ILS]
math: true
toc: true
---

## Introduction

- 단일 방정식 모형에서는 종속변수 $$Y$$와 독립변수 $$X$$ 간의 인과관계가 일방적($$X \rightarrow Y$$)이라고 가정한다. 그러나 경제 현상에서는 변수들이 서로 영향을 주고받는 양방향의 인과관계가 흔히 존재한다. 이처럼 두 개 이상의 endogenous variables가 여러 방정식을 통해 동시에 결정되는 체계를 **연립방정식 모형**이라고 한다.
- 연립방정식 체계 내의 개별 구조방정식(structural equation)을 OLS로 추정하면 Simultaneous Equation Bias가 발생한다. 이는 내생변수가 오차항과 상관관계를 가지기 때문이며, 결과적으로 OLS 추정량은 biased하며 inconsistent 하다.

## Identification Problem

- 연립방정식 모형을 적절히 추정하기 위해서는 먼저 구조방정식을 내생변수에 대해 푼 **축약형 방정식(Reduced-form)**으로 변환해야 한다. 축약형 방정식은 내생변수를 오직 exogenous variables와 오차항만의 함수로 표현한 것이다. 이 축약형 모형은 오차항과 독립변수 간의 상관관계가 없으므로 OLS로 일치성 있는 추정이 가능하다.
- **식별 문제**란, “관측 가능한 데이터로 추정한 reduced-form parameters로부터, 원래 우리가 알고 싶어 하는 structural parameters를 uniquely 역산해낼 수 있는가?”에 관한 것이다.

### 스칼라(Scalar) 버전

가장 고전적인 시장의 수요-공급 모형을 통해 스칼라 형태의 식별 문제를 확인해 보자.

**[구조 모형: 식별 불가능한 경우]**

$$
\begin{aligned}
\text{수요 방정식: } Q_t &= \alpha_0 + \alpha_1 P_t + u_{1t} \\
\text{공급 방정식: } Q_t &= \beta_0 + \beta_1 P_t + u_{2t}
\end{aligned}
$$

위 두 식에서 내생변수는 균형 거래량 $$Q_t$$와 균형 가격 $$P_t$$이다. 이를 $$P_t$$와 $$Q_t$$에 대해 풀면 다음과 같은 축약형 방정식을 얻는다.

$$
\begin{aligned}
P_t &= \frac{\beta_0 - \alpha_0}{\alpha_1 - \beta_1} + \frac{u_{2t} - u_{1t}}{\alpha_1 - \beta_1} \equiv \pi_0 + v_{1t} \\
Q_t &= \frac{\alpha_1 \beta_0 - \alpha_0 \beta_1}{\alpha_1 - \beta_1} + \frac{\alpha_1 u_{2t} - \beta_1 u_{1t}}{\alpha_1 - \beta_1} \equiv \pi_1 + v_{2t}
\end{aligned}
$$

우리가 관측 데이터를 통해 OLS로 안전하게 추정할 수 있는 값은 축약형 계수인 $$\pi_0$$와 $$\pi_1$$ 단 2개뿐이다. 하지만 우리가 찾아야 하는 구조 계수는 $$\alpha_0, \alpha_1, \beta_0, \beta_1$$ 등 4개이다. 2개의 방정식으로 4개의 미지수를 풀 수 없으므로, 이 모형은 **과소식별(Underidentified)** 상태이다.

**[외생변수 추가를 통한 식별]**
수요 방정식에 외생변수인 소득($$I_t$$)을 추가해 보자.

$$
\begin{aligned}
\text{수요 방정식: } Q_t &= \alpha_0 + \alpha_1 P_t + \alpha_2 I_t + u_{1t} \\
\text{공급 방정식: } Q_t &= \beta_0 + \beta_1 P_t + u_{2t}
\end{aligned}
$$

이 경우 축약형 방정식은 다음과 같이 외생변수 $$I_t$$를 포함하게 된다.

$$
\begin{aligned}
P_t &= \pi_0 + \pi_1 I_t + v_{1t} \\
Q_t &= \pi_2 + \pi_3 I_t + v_{2t}
\end{aligned}
$$

여기서 계수들의 관계를 정리하면 $$\pi_1 = \frac{-\alpha_2}{\alpha_1 - \beta_1}$$ 이고, $$\pi_3 = \frac{-\alpha_2 \beta_1}{\alpha_1 - \beta_1}$$ 이 된다. 흥미롭게도 두 계수를 나누면 다음이 성립한다.

$$
\frac{\pi_3}{\pi_1} = \beta_1
$$

즉, 수요 방정식에 포함된 외생변수 $$I_t$$ 덕분에 수요 곡선이 이동하게 되고, 이를 통해 고정된 공급 곡선의 궤적을 추적할 수 있게 된다. 그 결과 공급 방정식의 기울기 $$\beta_1$$을 **정확히 식별(Exactly identified)**할 수 있다. (단, 수요 방정식의 계수들은 여전히 식별되지 않는다.)

### Matrix 버전

연립방정식 모형을 $$M$$개의 내생변수와 $$K$$개의 외생(또는 선결)변수를 가진 일반적인 행렬 형태로 확장하면 다음과 같다.

$$
B Y_t + \Gamma X_t = U_t
$$

- $$Y_t$$: $$M \times 1$$ 차원의 내생변수 벡터, $$X_t$$: $$K \times 1$$ 차원의 외생변수 벡터
- $$U_t$$: $$M \times 1$$ 차원의 구조적 오차항 벡터, $$B$$: $$M \times M$$ 차원의 내생변수 계수 행렬
- $$\Gamma$$: $$M \times K$$ 차원의 외생변수 계수 행렬

&#8203;$$B$$가 non-singular matrix이라고 가정하고 양변의 앞에 $$B^{-1}$$을 곱하여 Reduced-form으로 변환하면 다음과 같다.

$$
\begin{aligned}
Y_t &= -B^{-1}\Gamma X_t + B^{-1}U_t \\
Y_t &= \Pi X_t + V_t
\end{aligned}
$$

여기서 축약형 계수 행렬 $$\Pi$$는 $$M \times K$$ 차원을 가지며, $$\Pi = -B^{-1}\Gamma$$ 가 성립한다.

우리는 OLS를 통해 $$\Pi$$ 행렬의 원소 $$M \times K$$개를 일치성 있게 추정할 수 있다. 그러나 알아내야 하는 구조 매개변수는 $$B$$ 행렬($$M^2$$개)과 $$\Gamma$$ 행렬($$MK$$개)에서 각 방정식의 정규화 조건(종속변수의 계수를 1로 설정하여 $$M$$개 차감)을 제외한 $$M^2 + MK - M$$개이다.

미지수의 수가 아는 값의 수보다 많으므로, 특정 변수를 방정식에서 제외하여 계수를 0으로 만드는 등의 사전적 제약(a priori restrictions) 없이는 $$B$$와 $$\Gamma$$를 도출할 수 없다.

이를 해결하기 위해 특정 $$i$$번째 방정식에 대해 다음의 두 가지 식별 조건이 충족되어야 한다.

#### 차수 조건 (Order Condition - 필요조건)

어떤 방정식이 식별되기 위해서는, 전체 모형에는 존재하지만 해당 방정식에서는 제외된 외생변수의 수($$K - k_i$$)가 해당 방정식에 포함된 내생변수의 수($$m_i$$)에서 1을 뺀 것보다 크거나 같아야 한다.

$$
K - k_i \ge m_i - 1
$$

- 적도식별 (Exactly identified)

  $$
  K - k_i = m_i - 1 \iff \text{Rank}(A) = M - 1
  $$

- 과잉식별 (Overidentified)

  $$
  K - k_i > m_i - 1 \iff \text{Rank}(A) = M - 1 \quad \& \quad \text{(배제된 변수의 수)} > M - 1
  $$

- 과소식별 (Underidentified, 식별 불가)

  $$
  K - k_i < m_i - 1 \iff \text{Rank}(A) < M
  $$

#### 계수 조건 (Rank Condition - 필요충분조건)

&#8203;$$i$$번째 방정식이 식별되기 위한 충분조건은, 전체 구조 계수 행렬 $$\begin{bmatrix} B & \Gamma \end{bmatrix}$$에서 해당 방정식에 **포함되지 않은(제외된)** 변수들에 해당하는 열들만 추출하여 만든 부분행렬의 Rank가 정확히 $$M - 1$$이 되어야 한다는 것이다.
이 두 가지 조건을 만족해야 ILS이나 2SLS 등을 활용하여 구조 계수를 성공적으로 추정할 수 있다.

## 추정법 간의 효율성 비교

설명의 편의를 위해 모든 변수는 평균이 0인 편차 형태로 변환되어 있다고 가정하며, 다음과 같은 단순한 구조방정식을 설정한다.

$$
y = \beta x + u
$$

여기서 $$y$$는 종속변수, $$x$$는 오차항 $$u$$와 상관관계를 가지는 내생 설명변수($$Cov(x, u) \neq 0$$)이다. 모형 외부에서 관측 가능한 외생변수(도구변수 후보)로 $$z_1, z_2$$ 등이 존재한다고 가정한다.

### Exactly Identified

도구변수를 $$z_1$$ 단 하나만 사용할 수 있는 상황이다. 내생변수가 1개이고 도구변수도 1개이므로 과대식별 상태이다. 이 상황에서는 세 가지 추정량이 대수적으로 완벽히 동일함을 증명할 수 있다.

#### ILS Estimator

ILS는 축약형 방정식을 먼저 추정한 뒤 구조 계수를 역산한다. $$y$$와 $$x$$를 도구변수 $$z_1$$에 대해 표현한 축약형 방정식은 다음과 같다.

$$
\begin{aligned}
x &= \pi_1 z_1 + v_1 \\
y &= \pi_2 z_1 + v_2
\end{aligned}
$$

각 식을 OLS로 추정하면 $$\hat{\pi}_1 = \frac{\sum z_1 x}{\sum z_1^2}$$ 이고 $$\hat{\pi}_2 = \frac{\sum z_1 y}{\sum z_1^2}$$ 이다.
ILS 추정량 $$\hat{\beta}_{ILS}$$는 두 축약형 계수의 비율로 도출된다.

$$
\hat{\beta}_{ILS} = \frac{\hat{\pi}_2}{\hat{\pi}_1} = \frac{\frac{\sum z_1 y}{\sum z_1^2}}{\frac{\sum z_1 x}{\sum z_1^2}} = \frac{\sum z_1 y}{\sum z_1 x}
$$

#### IV Estimator

단순 IV Estimator은 오차항과 도구변수의 표본 공분산이 0이 된다는 적률 조건($$\sum z_1 u = 0$$)을 이용하여 다음과 같이 구한다.

$$
\hat{\beta}_{IV} = \frac{\sum z_1 y}{\sum z_1 x}
$$

이를 통해 $$\hat{\beta}_{ILS} = \hat{\beta}_{IV}$$ 임을 바로 알 수 있다.

#### 2SLS Estimator

2SLS는 2단계에 걸쳐 추정한다. 먼저 1단계에서 내생변수 $$x$$를 도구변수 $$z_1$$에 대해 OLS 회귀한다. 이때 OLS 예측값은 $$\hat{x} = \hat{\pi}_1 z_1$$ 이 된다. (여기서 $$\hat{\pi}_1 = \frac{\sum z_1 x}{\sum z_1^2}$$) 그 다음 2단계에서 원래 방정식의 종속변수 $$y$$를 1단계의 예측값 $$\hat{x}$$에 대해 OLS 회귀한다.

$$
\hat{\beta}_{2SLS} = \frac{\sum \hat{x} y}{\sum \hat{x}^2}
$$

OLS 정규방정식의 직교성 성질에 의해 $$\sum \hat{x}^2 = \sum \hat{x} x$$ 가 성립하므로 식을 변형할 수 있다.

$$
\hat{\beta}_{2SLS} = \frac{\sum \hat{x} y}{\sum \hat{x} x}
$$

이제 $$\hat{x} = \hat{\pi}_1 z_1$$을 위 식에 대입하여 전개한다.

$$
\hat{\beta}_{2SLS} = \frac{\sum (\hat{\pi}_1 z_1) y}{\sum (\hat{\pi}_1 z_1) x} = \frac{\hat{\pi}_1 \sum z_1 y}{\hat{\pi}_1 \sum z_1 x}
$$

분자와 분모에서 공통 상수인 $$\hat{\pi}_1$$을 밖으로 빼내어 약분하면 다음과 같다.

$$
\hat{\beta}_{2SLS} = \frac{\sum z_1 y}{\sum z_1 x}
$$

적도 식별에서는 $$\hat{\beta}_{ILS} = \hat{\beta}_{IV} = \hat{\beta}_{2SLS}$$ 가 성립한다. 즉, 세 방법은 효율성 측면에서 아무런 차이가 없다.

### Overidentified

가용 도구변수가 $$z_1, z_2$$ 두 개 이상으로, 내생변수의 수보다 많아진 Overidentified 상황이다. 이 상황에서 ILS는 어떤 도구변수를 기준으로 식을 푸느냐에 따라 $$\hat{\beta}$$에 대한 해가 여러 개 도출되어 고유한 추정량을 구할 수 없으므로(추정 불가), 효율성 비교는 IV와 2SLS를 중심으로 이루어진다.

#### 단순 IV Estimator의 분산 도출

어떤 임의의 단일 도구변수 $$z$$ (예: $$z_1$$ 단독, 또는 임의로 가중평균한 $$0.5z_1 + 0.5z_2$$ 등)를 사용할 때, $$\hat{\beta}_{IV} = \frac{\sum z y}{\sum z x}$$ 의 점근적 분산은 오차항의 동분산성을 가정할 때 다음과 같이 근사된다.

$$
Var(\hat{\beta}_{IV}) \approx \frac{\sigma^2 \sum z^2}{(\sum z x)^2}
$$

분모와 분자를 표본의 크기 $$n$$으로 적절히 나누어 모분산과 모공분산의 형태로 표현하면 다음과 같다.

$$
Var(\hat{\beta}_{IV}) = \frac{\sigma^2}{n} \frac{Var(z)}{[Cov(z, x)]^2}
$$

여기서 두 변수의 상관계수 정의를 활용하여 수식을 변형해 보자.

$$
\begin{aligned}
Var(\hat{\beta}_{IV}) &= \frac{\sigma^2}{n} \frac{Var(z)}{[\rho_{zx} \sqrt{Var(z)Var(x)}]^2} \quad (\because Cov(z, x) = \rho_{zx} \sqrt{Var(z)Var(x)}) \\
&= \frac{\sigma^2}{n} \frac{Var(z)}{\rho_{zx}^2 Var(z) Var(x)} \\
&= \frac{\sigma^2}{n \cdot Var(x)} \cdot \frac{1}{\rho_{zx}^2}
\end{aligned}
$$

이 식은 매우 중요한 직관적 결론을 제시한다. 오차 분산 $$\sigma^2$$과 $$x$$의 분산 $$Var(x)$$가 주어진 상태에서, **추정량의 분산은 오직 도구변수 $$z$$와 내생변수 $$x$$ 간의 상관계수 제곱($$\rho_{zx}^2$$)에 반비례**한다. 즉, 도구변수가 내생변수를 잘 설명할수록(상관관계가 높을수록) 분산이 작아져 추정의 효율성이 높아진다.

#### 2SLS의 상관계수 극대화 원리

단순 IV는 넘치는 도구변수 중 일부만 사용하므로 상관계수 $$\rho_{zx}$$를 극대화하지 못하고 정보의 손실을 겪는다. 반면 2SLS는 1단계에서 $$x$$를 가용한 모든 도구변수 $$z_1, z_2$$에 대해 다중회귀분석한다.

$$
x = \pi_1 z_1 + \pi_2 z_2 + v
$$

이 1단계 회귀를 통해 예측값 $$\hat{x} = \hat{\pi}_1 z_1 + \hat{\pi}_2 z_2$$를 얻는다.
다중회귀분석의 OLS 원리상, $$\hat{x}$$는 관측된 $$z_1$$과 $$z_2$$의 선형결합으로 만들 수 있는 수많은 가상의 변수들 중에서 **실제 종속변수 $$x$$와의 상관계수를 극대화한 최적의 변수**이다.
즉, 임의로 조합한 어떤 단일 도구변수 $$z^*$$를 가져오더라도 다음의 부등식이 반드시 성립한다.

$$
\rho_{\hat{x} x}^2 \ge \rho_{z^* x}^2
$$

2SLS는 2단계에서 이 최적화된 $$\hat{x}$$를 단일 도구변수처럼 사용하는 것과 대수적으로 동일하다. 앞서 도출한 IV 분산 공식에 이를 대입하여 단순 IV와 비교하면 다음과 같다.

$$
Var(\hat{\beta}_{2SLS}) = \frac{\sigma^2}{n \cdot Var(x)} \cdot \frac{1}{\rho_{\hat{x} x}^2} \le \frac{\sigma^2}{n \cdot Var(x)} \cdot \frac{1}{\rho_{z^* x}^2} = Var(\hat{\beta}_{IV})
$$

결과적으로, $$Var(\hat{\beta}_{2SLS}) \le Var(\hat{\beta}_{IV})$$ 이 성립한다.

과잉식별 상태에서는 도구변수들의 가중치를 최적화하여 내생변수와의 상관관계를 극대화하는 **2SLS가 가용한 도구변수의 일부만 자의적으로 사용하는 단순 IV에 비해 추정량의 분산이 작아 점근적으로 우월한 효율성(Asymptotic Efficiency)**을 가진다.
