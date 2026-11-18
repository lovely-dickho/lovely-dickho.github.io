---
title: "[Econometrics] Limited Dependent Variable Models"
date: 2026-11-18 16:22:00 +0900
categories: [Econometrics, Notes]
tags: [LPM, Probit, Logit, Tobit, 다항로짓, 서열선다형]
math: true
toc: true
---

## 선형확률모형

선형확률모형은 종속변수 $$Y_i$$가 특정 사건의 발생 여부에 따라 1 또는 0의 값만을 가지는 binary 변수일 때, 이를 독립변수 $$X_i$$들의 선형결합으로 설명하는 회귀모형이다. 모형의 수식은 다음과 같이 정의된다.

$$
Y_i = \beta_1 + \beta_2 X_i + u_i
$$

&#8203;$$Y_i$$는 베르누이 확률분포를 따르는 확률변수이므로, 종속변수의 조건부 기댓값을 구하면 해당 사건이 발생할($$Y_i=1$$) 확률 그 자체가 된다. 이를 수식으로 나타내면 다음과 같다.

$$
\begin{aligned}
\mathbb{E}[Y_i | X_i] &= 1 \cdot P(Y_i = 1 | X_i) + 0 \cdot P(Y_i = 0 | X_i) \\
&= P(Y_i = 1 | X_i) \\
&= P_i
\end{aligned}
$$

즉, OLS를 이용해 추정한 조건부 기댓값 $$\mathbb{E}[Y_i | X_i]$$ 자체를 특정 사건이 발생할 조건부 확률 $$P_i$$로 해석하게 된다.

$$
\mathbb{E}[Y_i | X_i] = \beta_1 + \beta_2 X_i = P_i
$$

### 선형확률모형의 수리적, 통계적 한계점

LPM은 해석이 직관적이고 OLS로 간단히 추정할 수 있다는 장점이 있으나, 모형의 구조적 특성으로 인해 다음과 같은 계량경제학적 문제점들을 내포하고 있다.
#### 오차항의 비정규성

전통적인 고전적 선형회귀모형에서는 가설검정을 위해 오차항의 정규분포를 가정하지만, LPM의 오차항 $$u_i = Y_i - \beta_1 - \beta_2 X_i$$는 종속변수와 마찬가지로 이산적인 두 가지 확률값만을 가진다.

- $$Y_i = 1$$ 일 때: $$u_i = 1 - \beta_1 - \beta_2 X_i$$ (발생 확률: $$P_i$$)
- $$Y_i = 0$$ 일 때: $$u_i = - \beta_1 - \beta_2 X_i$$ (발생 확률: $$1-P_i$$)

따라서 오차항 $$u_i$$는 정규분포가 아닌 베르누이 분포를 따르며, 소표본에서의 통계적 검정(t-test, F-test 등)이 엄밀하게 성립하지 않는다. 다만 표본 크기가 커지면 중심극한정리(CLT)에 따라 OLS 추정량이 점근적 정규분포를 따르게 되므로 대표본 근사는 가능하다.

#### 오차항의 이분산성

&#8203;$$\mathbb{E}[u_i | X_i] = 0$$이라는 전제 하에 베르누이 분포를 따르는 오차항의 조건부 분산을 계산하면 다음과 같다.

$$
\begin{aligned}
Var(u_i | X_i) &= \mathbb{E}[u_i^2 | X_i] \\
&= (1 - \beta_1 - \beta_2 X_i)^2 P_i + (- \beta_1 - \beta_2 X_i)^2 (1 - P_i) \\
&= (1 - P_i)^2 P_i + (-P_i)^2 (1 - P_i) \\
&= P_i (1 - P_i)^2 + P_i^2 (1 - P_i) \\
&= P_i (1 - P_i) (1 - P_i + P_i) \\
&= P_i (1 - P_i)
\end{aligned}
$$

여기서 $$P_i = \beta_1 + \beta_2 X_i$$를 대입하면, 오차항의 분산은 아래와 같다.

$$
Var(u_i | X_i) = (\beta_1 + \beta_2 X_i)(1 - \beta_1 - \beta_2 X_i)
$$

즉, 오차항의 분산이 설명변수 $$X_i$$에 의존하는 구조가 된다. 이로 인해 OLS 추정량은 unbiasedness efficiency를 잃게 된다. 이를 보완하기 위해 $$w_i = \hat{P}_i(1-\hat{P}_i)$$를 가중치로 사용하는 WLS나 White test를 진행해볼 수 있다.

#### 조건부 확률의 논리적 범위 위배 (Nonfulfillment of 0 ≤ E[Yᵢ | Xᵢ] ≤ 1)

#### R²의 적합도 측도로서의 부적절성

실제 $$Y_i$$ 값들은 모두 0 또는 1의 y축 선상에만 존재하므로, 연속된 선형 회귀식으로 이 산포도를 완벽하게 관통하여 적합하는 것은 불가능에 가깝다. 따라서 모형의 실제 예측력과 관계없이 통상적으로 계산된 $$R^2$$는 대개 0.2$$\sim$$0.6 정도로 매우 낮게 도출된다

## 프로비트 모형과 로짓 모형

종속변수 $$Y_i$$가 0 또는 1의 값만 갖는 경우를 다룬다. 관측되지 않는 잠재변수 $$Y_i^* = \beta_0 + \beta_1 X_i + u_i$$를 가정하며, $$Y_i^* > 0$$일 때 $$Y_i = 1$$이 관측된다고 본다.

$$
P_i = P(Y_i=1|X_i) = P(u_i > -(\beta_0 + \beta_1 X_i)) = F(\beta_0 + \beta_1 X_i)
$$

오차항 $$u_i$$의 누적분포함수(CDF)인 $$F(\cdot)$$의 가정에 따라 두 모형으로 나뉜다.

#### 로짓 모형

&#8203;$$u_i$$가 표준 로지스틱 분포를 따른다.

$$
P_i = \Lambda(\beta_0 + \beta_1 X_i) = \frac{e^{\beta_0 + \beta_1 X_i}}{1 + e^{\beta_0 + \beta_1 X_i}}
$$

#### 프로빗 모형

&#8203;$$u_i$$가 표준 정규분포를 따른다.

$$
P_i = \Phi(\beta_0 + \beta_1 X_i) = \int_{-\infty}^{\beta_0 + \beta_1 X_i} \frac{1}{\sqrt{2\pi}} e^{-\frac{z^2}{2}} dz
$$

### MLE 도출

각 관측치는 0 또는 1을 갖는 독립적인 베르누이 시행이므로, likelihood function $$L$$과 $$\ln L$$은 다음과 같다.

$$
L(\beta_0, \beta_1) = \prod_{i=1}^n P_i^{Y_i} (1-P_i)^{1-Y_i}
$$

$$
\ln L = \sum_{i=1}^n \left[ Y_i \ln P_i + (1-Y_i) \ln(1-P_i) \right]
$$

이 식에 각 모형의 $$P_i$$를 대입한 후, $$\beta_0, \beta_1$$에 대해 편미분하여 1계 조건을 만족하는 파라미터를 수치적으로 추정한다.

### Marginal Effect 도출

확률 $$P_i$$는 $$X_i$$에 대해 비선형이므로, $$X_i$$가 확률에 미치는 한계효과는 연쇄법칙을 사용하여 도출한다. 잠재지수를 $$z = \beta_0 + \beta_1 X_i$$라 하면, $$P_i = F(z)$$이다.

$$
\frac{\partial P_i}{\partial X_i} = \frac{d F(z)}{d z} \cdot \frac{\partial z}{\partial X_i} = f(z) \cdot \beta_1
$$

여기서 $$f(z)$$는 해당 분포의 확률밀도함수(PDF)이다.

로짓 모형의 확률 정의는 $$P_i = \Lambda(z_i) = \frac{\exp(z_i)}{1+\exp(z_i)}$$ 이며, 여기서 $$z_i = \beta_0 + \beta_1 X_i$$이다.

$$
\frac{\partial P_i}{\partial X_i} = \frac{d \Lambda(z_i)}{d z_i} \cdot \frac{\partial z_i}{\partial X_i} \quad (\text{by Chain Rule})
$$

로지스틱 함수의 미분 특성 $$\frac{d \Lambda(z)}{dz} = \Lambda(z)(1-\Lambda(z))$$를 이용하면 다음과 같다.

$$
\frac{\partial P_i}{\partial X_i} = \Lambda(z_i)(1-\Lambda(z_i)) \cdot \beta_1 = P_i(1-P_i)\beta_1
$$

프로빗 모형의 확률 정의는 $$P_i = \Phi(z_i) = \int_{-\infty}^{z_i} \phi(t) dt$$ 이며, $$z_i = \beta_0 + \beta_1 X_i$$이다.

$$
\frac{\partial P_i}{\partial X_i} = \frac{d \Phi(z_i)}{d z_i} \cdot \frac{\partial z_i}{\partial X_i} \quad (\text{by Chain Rule})
$$

미적분학의 기본 정리에 의해 $$\frac{d \Phi(z)}{dz} = \phi(z)$$ (표준정규분포의 PDF)이므로 다음과 같다.

$$
\frac{\partial P_i}{\partial X_i} = \phi(\beta_0 + \beta_1 X_i) \cdot \beta_1
$$

두 모형 모두 한계효과는 **(해당 지점에서의 확률 밀도) $$\times$$ (회귀 계수)**의 형태를 띤다.
밀도함수 값($$P_i(1-P_i)$$ 또는 $$\phi(z_i)$$)은 항상 양수이므로, **한계효과의 부호는 오직 계수 $$\beta_1$$에 의해 결정**된다.

#### 계수 해석

계수 $$\beta_1$$ 자체는 한계효과가 아니다. $$\beta_1$$의 부호는 확률 증감의 방향을 나타내며, 실제 한계효과(크기)를 구하려면 반드시 표본 평균 등 특정 위치에서의 확률밀도($$f(z)$$)를 곱해주어야 한다.

## Tobit Model

종속변수가 0과 같은 하한선에서 절단(censored)되어 다수의 관측치가 0에 뭉쳐 있는 경우에 사용한다. 잠재변수 $$Y_i^* = \beta_0 + \beta_1 X_i + u_i$$ 이고, $$u_i \sim N(0, \sigma^2)$$를 가정한다.

$$
Y_i = \max(0, Y_i^*) = \begin{cases} \beta_0 + \beta_1 X_i + u_i & \text{if } Y_i^* > 0 \\ 0 & \text{if } Y_i^* \le 0 \end{cases}
$$

### MLE 도출

데이터가 0인 집단(censored)과 0보다 큰 양수인 집단(uncensored)으로 확률 및 밀도를 나누어 우도함수를 결합한다.

$$
L = \prod_{Y_i=0} \left[ 1 - \Phi\left(\frac{\beta_0 + \beta_1 X_i}{\sigma}\right) \right] \prod_{Y_i>0} \left[ \frac{1}{\sigma} \phi\left(\frac{Y_i - \beta_0 - \beta_1 X_i}{\sigma}\right) \right]
$$

이에 자연로그를 취하여 로그 우도함수를 구한다.

$$
\ln L = \sum_{Y_i=0} \ln \left[ 1 - \Phi\left(\frac{\beta_0 + \beta_1 X_i}{\sigma}\right) \right] + \sum_{Y_i>0} \left[ -\ln \sigma + \ln \phi\left(\frac{Y_i - \beta_0 - \beta_1 X_i}{\sigma}\right) \right]
$$

### 한계효과와 계수 해석

#### 잠재변수(Yᵢ∗)에 대한 효과

&#8203;$$\mathbb{E}[Y_i^*|X_i] = \beta_0 + \beta_1 X_i$$ 이므로, $$\beta_1$$은 관측되지 않는 잠재변수에 대한 한계효과를 의미한다.

$$
\frac{\partial \mathbb{E}[Y_i^*|X_i]}{\partial X_i} = \beta_1
$$

#### 실제 관측치(Yᵢ)에 대한 효과

실제 기댓값 $$\mathbb{E}[Y_i|X_i]$$는 밀스비(Mills ratio) 성질에 따라 $$\mathbb{E}[Y_i|X_i] = \Phi(z)\sigma z + \sigma\phi(z)$$ 로 구성되며(단, $$z = \frac{\beta_0 + \beta_1 X_i}{\sigma}$$), 이를 미분하면 연쇄법칙에 의해 복잡한 항이 소거되고 다음과 같이 축소된 형태가 도출된다.

$$
\frac{\partial \mathbb{E}[Y_i|X_i]}{\partial X_i} = \Phi\left(\frac{\beta_0 + \beta_1 X_i}{\sigma}\right) \cdot \beta_1
$$

즉, 토빗 모형에서 실제 관측 데이터에 미치는 한계효과는 잠재변수의 계수($$\beta_1$$)에 ‘절단되지 않고 관측될 확률(정규분포 CDF)’을 곱한 값으로 해석한다. 0에 데이터가 뭉쳐있기 때문에 한계효과가 $$\beta_1$$보다 항상 작게(축소되어) 나타난다.

## 다항 로지트 모형

다항 로지트 모형은 개인 $$i$$의 특성이 대안 선택에 미치는 영향을 분석한다. 개인의 특성을 나타내는 $$K$$개의 설명변수가 있다고 가정하자.

#### 효용 함수:

개인 $$i$$가 대안 $$j$$로부터 얻는 확정적 효용 $$V_{ij}$$는 다음과 같이 각 설명변수 $$X_{ki}$$와 대안별로 다른 계수 $$\beta_{kj}$$의 선형 결합으로 나타낸다. 상수항 $$\beta_{0j}$$도 대안별로 존재한다.

$$
V_{ij} = \beta_{0j} + \beta_{1j} X_{1i} + \beta_{2j} X_{2i} + \dots + \beta_{Kj} X_{Ki} = \beta_{0j} + \sum_{k=1}^K \beta_{kj} X_{ki}
$$

#### 선택 확률:

개인 $$i$$가 $$J$$개의 대안 중 대안 $$j$$를 선택할 확률 $$P_{ij}$$는 오차항이 제1종 극단값 분포를 따른다는 가정 하에 다음과 같이 도출된다.

$$
P_{ij} = \frac{\exp\left(\beta_{0j} + \sum_{k=1}^K \beta_{kj} X_{ki}\right)}{\sum_{m=1}^J \exp\left(\beta_{0m} + \sum_{k=1}^K \beta_{km} X_{ki}\right)}
$$

#### 모형의 식별 (정규화):

계수의 유일한 해를 구하기 위해 기저 대안(예: $$j=1$$)의 모든 파라미터를 0으로 제약한다 ($$\beta_{01} = 0$$, $$\beta_{k1} = 0 \text{ for } k=1,\dots,K$$).

$$
\begin{aligned}
P_{i1} &= \frac{1}{1 + \sum_{m=2}^J \exp\left(\beta_{0m} + \sum_{k=1}^K \beta_{km} X_{ki}\right)} \\
P_{ij} &= \frac{\exp\left(\beta_{0j} + \sum_{k=1}^K \beta_{kj} X_{ki}\right)}{1 + \sum_{m=2}^J \exp\left(\beta_{0m} + \sum_{k=1}^K \beta_{km} X_{ki}\right)} \quad (j \neq 1)
\end{aligned}
$$

## 조건부 로지트 모형

조건부 로지트 모형은 각 대안이 가진 고유한 특성이 선택에 미치는 영향을 분석한다. 대안의 속성을 나타내는 $$M$$개의 설명변수가 있다고 가정한다.

#### 효용 함수:

대안 $$j$$의 특성 변수 $$Z_{mij}$$에 대한 사람들의 한계효용(파라미터) $$\gamma_m$$은 대안의 종류에 관계없이 동일하다고 가정한다.

$$
V_{ij} = Z_{1ij} \gamma_1 + Z_{2ij} \gamma_2 + \dots + Z_{Mij} \gamma_M = \sum_{m=1}^M \gamma_m Z_{mij}
$$

(상수항은 모형의 확률식에서 분모와 분자에 공통으로 들어가 약분되므로 식별되지 않아 제외한다.)

#### 선택 확률:

$$
P_{ij} = \frac{\exp\left(\sum_{m=1}^M \gamma_m Z_{mij}\right)}{\sum_{l=1}^J \exp\left(\sum_{m=1}^M \gamma_m Z_{mil}\right)}
$$

## 서열선다형 모형

서열선다형 모형은 관측되지 않는 연속형 잠재변수 $$y_i^*$$를 도입하여 설명한다. 개인 $$i$$의 성향이나 효용을 나타내는 이 잠재변수를 $$K$$개의 개별 설명변수들의 선형 결합으로 가정한다.

$$
y_i^* = \beta_1 X_{1i} + \beta_2 X_{2i} + \dots + \beta_K X_{Ki} + \varepsilon_i = \sum_{k=1}^K \beta_k X_{ki} + \varepsilon_i
$$

여기서 $$X_{ki}$$는 개인 $$i$$의 $$k$$번째 설명변수, $$\beta_k$$는 그에 대응하는 파라미터, $$\varepsilon_i$$는 오차항이다. 통상적으로 식별을 위해 상수항 $$\beta_0$$는 제외하거나, 임계치 중 하나를 0으로 제약한다.

### 관측변수와 임계치

우리가 실제 데이터에서 관측하는 변수 $$y_i$$는 $$J$$개의 서열화된 범주($$1, 2, \dots, J$$) 중 하나이다. 잠재변수 $$y_i^*$$가 어느 구간에 속하느냐에 따라 관측되는 $$y_i$$가 결정되며, 이 구간을 나누는 미지의 파라미터들을 임계치(thresholds) $$\mu_1 < \mu_2 < \dots < \mu_{J-1}$$라고 한다.

$$
y_i =
\begin{cases}
1 & \text{if } -\infty < \sum_{k=1}^K \beta_k X_{ki} + \varepsilon_i \le \mu_1 \\
2 & \text{if } \mu_1 < \sum_{k=1}^K \beta_k X_{ki} + \varepsilon_i \le \mu_2 \\
\vdots & \vdots \\
J & \text{if } \mu_{J-1} < \sum_{k=1}^K \beta_k X_{ki} + \varepsilon_i < \infty
\end{cases}
$$

편의상 $$\mu_0 = -\infty$$, $$\mu_J = \infty$$로 정의한다.

### 선택 확률

개인 $$i$$가 범주 $$j$$를 선택할 확률은 잠재변수 $$y_i^*$$가 $$\mu_{j-1}$$과 $$\mu_j$$ 사이에 존재할 확률이다. 오차항 $$\varepsilon_i$$의 CDF를 $$F(\cdot)$$라고 할 때, 스칼라 형태의 확률은 다음과 같이 도출된다.

$$
\begin{aligned}
P(y_i = j) &= P\left(\mu_{j-1} < \sum_{k=1}^K \beta_k X_{ki} + \varepsilon_i \le \mu_j\right) \\
P(y_i = j) &= P\left(\mu_{j-1} - \sum_{k=1}^K \beta_k X_{ki} < \varepsilon_i \le \mu_j - \sum_{k=1}^K \beta_k X_{ki}\right) \\
P(y_i = j) &= F\left(\mu_j - \sum_{k=1}^K \beta_k X_{ki}\right) - F\left(\mu_{j-1} - \sum_{k=1}^K \beta_k X_{ki}\right)
\end{aligned}
$$

오차항이 표준 정규분포를 따르면 $$F(\cdot) = \Phi(\cdot)$$가 되어 서열 프로빗(Ordered Probit) 모형이 되고, 로지스틱 분포를 따르면 $$F(\cdot) = \Lambda(\cdot)$$가 되어 서열 로지트(Ordered Logit) 모형이 된다.

### 한계효과 및 해석의 주의점

서열선다형 모형에서 계수 $$\beta_m$$의 부호가 의미하는 바를 해석할 때는 매우 주의해야 한다. 특정 설명변수 $$X_{mi}$$가 1단위 증가할 때 특정 범주 $$j$$에 속할 확률의 변화인 한계효과는 위 확률식을 $$X_{mi}$$에 대해 편미분하여 구한다. $$F(\cdot)$$의 확률밀도함수(PDF)를 $$f(\cdot)$$라고 할 때, chain rule에 의해 한계효과는 다음과 같이 전개된다.

$$
\frac{\partial P(y_i = j)}{\partial X_{mi}} = \left[ f\left(\mu_{j-1} - \sum_{k=1}^K \beta_k X_{ki}\right) - f\left(\mu_j - \sum_{k=1}^K \beta_k X_{ki}\right) \right] \beta_m
$$

이 수식을 양 끝단 범주와 중간 범주로 나누어 살펴보면 경제학적 시사점이 명확해진다.

#### 1) 가장 낮은 범주 (j=1):

&#8203;$$\mu_0 = -\infty$$ 이므로 $$f(-\infty) = 0$$ 이다.

$$
\frac{\partial P(y_i = 1)}{\partial X_{mi}} = - f\left(\mu_1 - \sum_{k=1}^K \beta_k X_{ki}\right) \beta_m
$$

확률밀도함수 $$f(\cdot)$$는 항상 양수이므로, 만약 $$\beta_m > 0$$ 이라면 가장 낮은 범주에 속할 확률은 무조건 감소한다.

#### 2) 가장 높은 범주 (j=J):

&#8203;$$\mu_J = \infty$$ 이므로 $$f(\infty) = 0$$ 이다.

$$
\frac{\partial P(y_i = J)}{\partial X_{mi}} = f\left(\mu_{J-1} - \sum_{k=1}^K \beta_k X_{ki}\right) \beta_m
$$

마찬가지로 $$\beta_m > 0$$ 이라면 가장 높은 범주에 속할 확률은 무조건 증가한다.

#### 3) 중간 범주 (1 < j < J):

$$
\frac{\partial P(y_i = j)}{\partial X_{mi}} = \left[ f\left(\mu_{j-1} - \sum_{k=1}^K \beta_k X_{ki}\right) - f\left(\mu_j - \sum_{k=1}^K \beta_k X_{ki}\right) \right] \beta_m
$$

중간 범주의 경우, 대괄호 안의 두 확률밀도함수 값의 차이의 부호가 개인의 특성 $$X_i$$ 값과 임계치의 위치에 따라 달라진다. 따라서 $$\beta_m$$이 양수(+)라고 해서 중간 범주에 속할 확률이 반드시 증가하거나 감소한다고 단정 지을 수 없다.

결론적으로, 서열선다형 모형에서는 단순히 추정된 파라미터 $$\beta_m$$의 부호만으로 중간 대안들의 확률 변화 방향을 예측할 수 없으므로, 반드시 스칼라 편미분식을 통해 각 범주별 한계효과를 직접 계산하고 부호를 확인해야 한다.
