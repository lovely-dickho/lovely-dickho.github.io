---
title: "[Econometrics] Multicollinearity and Autocorrelation"
date: 2026-11-16 17:55:00 +0900
categories: [Econometrics]
tags: [다중공선성, 자기상관, Durbin-Watson, Breusch-Godfrey, Ljung-Box]
math: true
toc: true
---

## Multicolinearity

#### Scalar Version

다음의 모형을 가정해보자.

$$
y_i = \beta_0 + \beta_1 X_{1i} + \beta_2 X_{2i} + u_i
$$

이 모형에서 OLS를 통해 도출한 추정량 $$\hat{\beta}_1$$의 분산은 다음과 같이 계산된다.

$$
Var(\hat{\beta}_1) = \frac{\sigma^2}{\sum_{i=1}^{n} (X_{1i} - \bar{X}_1)^2 (1 - R_1^2)} = \frac{\sigma^2}{SST_1 (1 - R_1^2)}
$$

여기서 $$R_1^2$$는 보조 회귀분석의 결정계수이다. 즉, 독립변수 $$X_{1i}$$를 또 다른 독립변수 $$X_{2i}$$에 대해 회귀분석했을 때 도출되는 $$R^2$$ 값이다.

$$
X_{1i} = \gamma_0 + \gamma_1 X_{2i} + v_i \quad \implies \quad R_1^2
$$

- $$X_1$$과 $$X_2$$ 간에 상관관계가 전혀 없다면 단순회귀분석의 분산과 일치하게 된다.

  $$
  R_1^2 = 0  \longrightarrow Var(\hat{\beta}_1) = \frac{\sigma^2}{SST_1}
  $$

- $$X_1$$과 $$X_2$$ 간에 다중공선성가 존재한다면, $$R_1^2 \to 1$$에 가까워진다.
- $$R_1^2$$가 1에 가까워질수록 분모의 $$(1 - R_1^2)$$ 항이 0에 수렴하게 되므로, 분산 $$Var(\hat{\beta}_1)$$은 무한대로 폭발하게 된다. 이처럼 분산을 팽창시키는 항을 분산팽창지수(VIF)라 한다.

  $$
  \frac{1}{1 - R_1^2}
  $$

분산이 폭발적으로 커지면 *t*-통계량이 매우 작아져서, 실제로는 유의미한 변수임에도 불구하고 귀무가설($$H_0: \beta_1 = 0$$)을 기각하지 못하는 제2종 오류가 발생할 확률이 극도로 높아진다.

#### Matrix Version

&#8203;$$Y$$는 $$n \times 1$$, $$X$$는 $$n \times k$$, $$\beta$$는 $$k \times 1$$, $$U$$는 $$n \times 1$$ 행렬이다.

$$
Y = X\beta + U
$$

오차항의 동분산성 가정 하에서 OLS 추정량 $$\hat{\beta}$$의 분산-공분산 행렬은 다음과 같다.

$$
Var(\hat{\beta}) = \sigma^2 (X^\mathsf{T}X)^{-1}
$$

역행렬의 정의에 따라 위 식은 수반행렬과 행렬식을 이용해 다음과 같이 풀어서 쓸 수 있다.

$$
(X^\mathsf{T}X)^{-1} = \frac{\text{Adj}(X^\mathsf{T}X)}{\det(X^\mathsf{T}X)}
$$

#### Perfect Multicollinearity

독립변수들 사이에 완벽한 선형 관계가 존재하면, $$X$$ 행렬의 열벡터들이 선형 종속이 된다.
이 경우 행렬 $$X^\mathsf{T}X$$는 Non-singular가 아니게 되므로 행렬식 $$\det(X^\mathsf{T}X) = 0$$이 된다. 따라서 역행렬이 아예 존재하지 않아 추정량 $$\hat{\beta}$$ 자체를 계산할 수 없다.

#### High Multicollinearity

독립변수들이 완벽한 선형 관계는 아니지만 매우 강한 상관관계를 가질 때 발생한다. 이 경우 $$X$$ 행렬의 열벡터들이 ‘거의’ 선형 종속에 가깝다. 즉, 행렬식 $$\det(X^\mathsf{T}X)$$의 값이 0에 매우 가깝게 수렴($$\det(X^\mathsf{T}X) \approx 0$$)함을 의미한다. 결과적으로 분모가 0에 수렴하므로, 역행렬 $$(X^\mathsf{T}X)^{-1}$$의 대각 원소들이 기하급수적으로 커지게 된다.

## AutoCorrelation

### Durbin-Watson 검정

Durbin-Watson 검정은 회귀모형의 잔차를 이용하여 1계 자기상관의 존재 여부를 파악하는 가장 대표적인 방법이다. 오차항이 AR(1) 과정을 따른다고 가정한다.

$$
u_t = \rho u_{t-1} + \varepsilon_t
$$

단, 이 검정은 설명변수($$X$$)가 비확률적이며, 종속변수의 시차항이 설명변수에 포함되지 않았다는 가정 하에서만 유효하다.

**①** 주어진 모형을 OLS로 추정하여 잔차 $$\hat{u}_t$$를 구한다.

**②** 도출된 잔차를 바탕으로 다음과 같이 $$d$$ 검정 통계량을 정의한다.

$$
d = \frac{\sum_{t=2}^n (\hat{u}_t - \hat{u}_{t-1})^2}{\sum_{t=1}^n \hat{u}_t^2}
$$

위 식의 분자를 전개하면 다음과 같다.

$$
\sum_{t=2}^n (\hat{u}_t - \hat{u}_{t-1})^2 = \sum_{t=2}^n \hat{u}_t^2 + \sum_{t=2}^n \hat{u}_{t-1}^2 - 2 \sum_{t=2}^n \hat{u}_t \hat{u}_{t-1}
$$

**③** 표본의 크기 $$n$$이 충분히 크다면, $$\sum_{t=2}^n \hat{u}_t^2 \approx \sum_{t=2}^n \hat{u}_{t-1}^2 \approx \sum_{t=1}^n \hat{u}_t^2$$ 로 근사할 수 있다. 이를 $$d$$ 통계량 식에 대입하면 다음과 같은 관계가 도출된다.

$$
\begin{aligned}
d &\approx \frac{2 \sum_{t=1}^n \hat{u}_t^2 - 2 \sum_{t=2}^n \hat{u}_t \hat{u}_{t-1}}{\sum_{t=1}^n \hat{u}_t^2} \\
&= 2 \left( 1 - \frac{\sum_{t=2}^n \hat{u}_t \hat{u}_{t-1}}{\sum_{t=1}^n \hat{u}_t^2} \right)
\end{aligned}
$$

**④** 괄호 안의 분수 항은 표본 자기상관 계수 $$\hat{\rho}$$의 추정치와 동일하므로, 다음이 성립한다.

$$
d \approx 2(1 - \hat{\rho})
$$

**⑤** $$\hat{\rho}$$의 범위가 $$[-1, 1]$$이므로, $$d$$ 통계량의 범위는 $$[0, 4]$$가 된다. 이를 바탕으로 가설을 검정한다.

- $$\hat{\rho} = 0$$ (자기상관 없음) $$\implies d \approx 2$$
- $$\hat{\rho} \to 1$$ (강한 양(+)의 자기상관) $$\implies d \to 0$$
- $$\hat{\rho} \to -1$$ (강한 음(-)의 자기상관) $$\implies d \to 4$$

### Durbin h 검정

**회귀모형에 종속변수의 시차항이 설명변수로 포함된 경우**, OLS 추정량은 일치성을 잃으며 Durbin-Watson $$d$$ 통계량은 점근적으로 2에 편의되는 문제가 발생한다. Durbin $$h$$ 검정은 이를 극복하기 위해 고안된 대안적 검정 방법이다.
다음과 같은 자기회귀 모형을 가정한다.

$$
Y_t = \beta_1 + \beta_2 X_t + \gamma Y_{t-1} + u_t
$$

**①** 기본 OLS 회귀를 수행하여 $$d$$ 통계량을 구하고, 이를 바탕으로 표본 자기상관 계수 $$\hat{\rho} \approx 1 - d/2$$ 를 계산한다.

**②** 시차항 종속변수의 계수 추정량인 $$\hat{\gamma}$$의 분산 $$\widehat{\text{var}}(\hat{\gamma})$$를 추출한다.

**③** 도출된 값들을 이용하여 다음과 같이 Durbin $$h$$ 통계량을 정의하고 계산한다.

$$
h = \hat{\rho} \sqrt{\frac{n}{(1 - n) \cdot \widehat{\text{var}}(\hat{\gamma})}} \quad \sim_{asy} \quad  N(0, 1)
$$

여기서 $$n$$은 표본의 크기이다.

**④** 표본의 크기 $$n$$이 클 때, 귀무가설($$\rho = 0$$) 하에서 $$h$$ 통계량은 점근적으로 표준정규분포를 따른다.

**⑤** 유의수준 $$\alpha$$에서 계산된 $$|h|$$ 값이 임계치보다 크면 귀무가설을 기각하고 1계 자기상관이 존재한다고 판단한다.

### Breusch-Godfrey (BG) 검정 (LM Test)

BG 검정은 1계 자기상관만을 가정하는 Durbin-Watson 검정의 한계를 극복하여 고차(high-order) 자기상관을 검정할 수 있으며, 설명변수에 종속변수의 시차항이 포함되어 있어도 유효하게 사용할 수 있는 매우 일반적인 방법이다.
오차항 $$u_t$$가 $$AR(p)$$ 과정을 따른다고 가정한다.

$$
u_t = \rho_1 u_{t-1} + \rho_2 u_{t-2} + \cdots + \rho_p u_{t-p} + \varepsilon_t
$$

귀무가설($$H_0$$)은 어떠한 차수에서도 자기상관이 존재하지 않는다($$\rho_1 = \rho_2 = \cdots = \rho_p = 0$$)는 것이다.

**①** 원래의 회귀모형($$Y_t = \beta_1 + \beta_2 X_t + \cdots + u_t$$)을 OLS로 추정하여 잔차 $$\hat{u}_t$$를 추출한다.

**②** 이 잔차를 종속변수로 하고, 원래 모형의 모든 설명변수들과 잔차의 시차항들($$\hat{u}_{t-1}, \hat{u}_{t-2}, \cdots, \hat{u}_{t-p}$$)을 설명변수로 포함하는 다음의 보조 회귀식을 추정한다.

$$
\hat{u}_t = \alpha_1 + \alpha_2 X_t + \cdots + \hat{\rho}_1 \hat{u}_{t-1} + \hat{\rho}_2 \hat{u}_{t-2} + \cdots + \hat{\rho}_p \hat{u}_{t-p} + v_t
$$

이 보조 회귀식으로부터 결정계수 $$R^2$$를 도출한다.

**③** 귀무가설이 참일 때, 표본 크기가 충분히 크다면 다음의 검정 통계량은 점근적으로 자유도가 $$p$$인 카이제곱($$\chi^2$$) 분포를 따른다.

$$
(n-p)R^2 \sim_{asy} \chi^2_p
$$

**④** 계산된 $$(n-p)R^2$$ 통계량이 유의수준의 카이제곱 임계치를 초과하면 귀무가설을 기각하고 적어도 하나 이상의 고차 자기상관이 존재한다고 결론 내린다.

### Ljung-Box (LB) 검정

Ljung-Box 검정은 데이터 자체나 회귀모형의 잔차가 white noise인지, 즉 특정 시차 $$m$$까지의 자기상관 계수들이 결합적으로 0인지 검정할 때 사용된다. 귀무가설($$H_0$$)은 시차 $$1$$부터 $$m$$까지의 모든 자기상관 계수가 0이라는 것이다 ($$\rho_1 = \rho_2 = \cdots = \rho_m = 0$$).

**①** 분석하고자 하는 잔차 $$\hat{u}_t$$의 시차 $$k$$에 대한 표본 자기상관 함수(SACF) $$\hat{\rho}_k$$를 구한다.

$$
\hat{\rho}_k = \frac{\sum_{t=k+1}^n \hat{u}_t \hat{u}_{t-k}}{\sum_{t=1}^n \hat{u}_t^2}
$$

**②** 표본 크기 $$n$$과 검정하고자 하는 최대 시차 $$m$$을 사용하여, Ljung-Box 통계량 $$LB$$를 다음과 같이 계산한다.

$$
LB = n(n+2) \sum_{k=1}^m \frac{\hat{\rho}_k^2}{n-k} \quad  \sim_{asy} \quad N(0, 1)
$$

여기서 분모의 $$n-k$$는 각 시차 $$k$$가 증가할수록 계산에 사용되는 데이터 쌍의 수가 줄어드는 것을 보정하는 역할을 한다.

**③** 귀무가설 하에서, $$LB$$ 통계량은 점근적으로 자유도가 $$m$$인 카이제곱 분포를 따른다. 단, ARMA($$p, q$$) 모형의 잔차에 대해 검정할 경우, 자유도는 $$m-p-q$$가 된다.

**④** 도출된 $$LB$$ 통계량이 임계치보다 크면 귀무가설을 기각하여 시차 $$m$$ 이내에 자기상관이 존재한다고 판단한다.
