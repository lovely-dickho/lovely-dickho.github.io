---
title: "[Econometrics] Omitted Variable Problem"
date: 2026-11-16 17:56:00 +0900
categories: [Econometrics, Notes]
tags: [OVP, 누락변수, 측정오차, Hausman, RESET, Chow]
math: true
toc: true
---

## Omitted Variable Problem

분석에 반드시 포함되어야 할 중요한 변수가 모형에서 누락되었을 때, 이 누락된 변수가 포함된 독립변수와 상관관계를 가진다면 추정량에 Bias가 발생한다.

True Model이 다음과 같다고 가정하자.

$$
y_i = \beta_0 + \beta_1 x_{1i} + \beta_2 x_{2i} + u_i
$$

하지만 $$x_{2i}$$를 관측할 수 없어 누락하고, 다음과 같은 모형을 OLS로 추정한다고 하자.

$$
y_i = \tilde{\beta}_0 + \tilde{\beta}_1 x_{1i} + \tilde{u}_i \quad (\text{단, } \tilde{u}_i = \beta_2 x_{2i} + u_i)
$$

이때 $$\tilde{\beta}_1$$의 OLS 추정량 식에 $$y_i$$를 대입하여 기댓값을 취하면 다음과 같이 편의가 도출된다.

$$
\begin{aligned}
\tilde{\beta}_1 &= \frac{Cov(x_1, y)}{Var(x_1)} = \beta_1 + \beta_2 \frac{Cov(x_1, x_2)}{Var(x_1)} + \frac{Cov(x_1, u)}{Var(x_1)} \\
\mathbb{E}[\tilde{\beta}_1] &= \beta_1 + \beta_2 \frac{Cov(x_1, x_2)}{Var(x_1)} \quad (\because \mathbb{E}[Cov(x_1, u)] = 0)
\end{aligned}
$$

여기서 $$\tilde{\delta}_1 = \frac{Cov(x_1, x_2)}{Var(x_1)}$$는 $$x_2$$를 $$x_1$$에 대해 회귀분석했을 때의 기울기이다. 즉, 누락변수 편의의 크기는 누락된 변수의 진정한 효과($$\beta_2$$)와 포함된 변수와 누락된 변수 간의 상관관계($$\tilde{\delta}_1$$)의 곱으로 결정된다.

$$
Bias(\tilde{\beta}_1) = \mathbb{E}[\tilde{\beta}_1] - \beta_1 = \beta_2 \tilde{\delta}_1
$$

## Detection for Endogeneity

### 램지 RESET 검정

RESET 검정은 누락된 변수나 잘못된 함수 형태 등으로 인해 모형에 전반적인 specification error가 존재하는지 확인하는 검정이다.

- **귀무가설($$H_0$$):** $$\gamma_1 = \gamma_2 = \cdots = \gamma_m = 0$$ (모형의 설정 오류가 없다. 즉, 추가된 비선형 항들이 유의하지 않다.)
- **대립가설($$H_1$$):** 적어도 하나의 $$\gamma_j \neq 0$$ (모형에 설정 오류가 존재한다.)

**검정 절차:**

**①** 원래의 기본 회귀모형(제약 모형)을 OLS로 추정하여 제한된 잔차제곱합 $$RSS_R$$을 구한다.

$$
Y_i = \beta_1 + \beta_2 X_{2i} + \cdots + \beta_k X_{ki} + u_i
$$

**②** 이 모형으로부터 구한 종속변수의 적합값 $$\hat{Y}_i$$의 거듭제곱 항들($$\hat{Y}_i^2, \hat{Y}_i^3$$ 등 $$m$$개)을 새로운 설명변수로 추가하여 보조 회귀식(무제약 모형)을 구성하고, OLS로 추정하여 무제약 잔차제곱합 $$RSS_{UR}$$을 구한다.

$$
Y_i = \beta_1 + \beta_2 X_{2i} + \cdots + \beta_k X_{ki} + \gamma_1 \hat{Y}_i^2 + \gamma_2 \hat{Y}_i^3 + \cdots + \gamma_m \hat{Y}_i^{m+1} + v_i
$$

**③** 두 모형의 잔차제곱합을 이용하여 다음의 $$F$$ 검정 통계량을 계산한다.

$$
F = \frac{(RSS_R - RSS_{UR}) / m}{RSS_{UR} / (n - k - m)}
$$

여기서 $$n$$은 표본 크기, $$k$$는 원래 모형의 파라미터(상수항 포함) 개수, $$m$$은 새롭게 추가된 거듭제곱 항의 개수이다.

**④** 계산된 $$F$$ 통계량이 유의수준의 $$F_{(m, n-k-m)}$$ 임계치를 초과하면 귀무가설을 기각하여 설정 오류가 존재한다고 판단한다.

### Chow 검정

Chow 검정은 특정 시점이나 기준을 바탕으로 표본을 나누었을 때, coefficient에 구조적 변화가 발생했는지 검정한다.

- **귀무가설($$H_0$$):** $$\boldsymbol{\beta}_1 = \boldsymbol{\beta}_2$$ (구조적 변화가 없다)
- **대립가설($$H_1$$):** $$\boldsymbol{\beta}_1 \neq \boldsymbol{\beta}_2$$ (구조적 변화가 존재한다.)

**검정 절차:**
전체 표본 크기가 $$n = n_1 + n_2$$이고, 추정할 모수의 개수가 $$k$$라고 가정한다.

**①** 전체 표본 $$n$$개로 pooled OLS 회귀를 수행하여 제한된 잔차제곱합 $$RSS_R$$을 구한다.

**②** 표본을 $$n_1$$과 $$n_2$$ 두 그룹으로 분리하여 각각 OLS를 수행하고, 첫 번째 그룹의 잔차제곱합 $$RSS_1$$과 두 번째 그룹의 잔차제곱합 $$RSS_2$$를 더해 무제약 잔차제곱합 $$URSS = RSS_1 + RSS_2$$를 구한다.

**③** 다음의 $$F$$ 통계량을 통해 가설을 검정한다.

$$
F = \frac{(RSS_R - URSS) / k}{URSS / (n_1 + n_2 - 2k)} \sim F_{(k, n_1 + n_2 - 2k)}
$$

**④** $$F$$ 통계량이 임계치를 초과하면 귀무가설을 기각하여 구조적 변화가 존재한다고 결론 내린다.

### 하우스만 검정 (Hausman Test / Durbin-Wu-Hausman Test)

설명변수와 오차항 간의 상관관계, 즉 endogeneity가 존재하는지 확인한다.

- **귀무가설($$H_0$$):** $$Cov(X, u) = 0$$ (내생성이 존재하지 않는다. OLS 추정량이 일치성과 효율성을 갖는다.)
- **대립가설($$H_1$$):** $$Cov(X, u) \neq 0$$ (내생성이 존재한다. IV 추정량만이 일치성을 갖는다.)

우선 스칼라 버전으로 접근해보자. 관심 있는 구조방정식이 $$Y_i = \beta_0 + \beta_1 X_i + \beta_2 W_i + u_i$$ 이고, 내생 의심 변수 $$X_i$$에 대한 도구변수 $$Z_i$$가 존재할 때 적용한다.

**①** 내생 의심 변수 $$X_i$$를 모든 외생변수 $$W_i$$와 도구변수 $$Z_i$$에 대해 회귀한다.

$$
X_i = \alpha_0 + \alpha_1 Z_i + \alpha_2 W_i + v_i
$$

이 회귀식으로부터 잔차 $$\hat{v}_i$$를 추출한다.

**②** 원래의 구조방정식에 1단계에서 구한 잔차 $$\hat{v}_i$$를 설명변수로 추가하여 OLS를 수행한다.

$$
Y_i = \beta_0 + \beta_1 X_i + \beta_2 W_i + \delta \hat{v}_i + e_i
$$

**③** 추가된 잔차항의 계수 $$\delta$$에 대해 통상적인 $$t$$ 검정을 수행한다 ($$H_0: \delta = 0$$). $$\delta$$가 통계적으로 유의미하면 귀무가설을 기각하고 $$X_i$$에 내생성이 존재한다고 결론 내린다.

이제 행렬로 접근해보자. 다중 내생변수 및 파라미터 벡터 전체에 대해 일반화된 형태로 접근한다.

**①** 모수 벡터 $$\boldsymbol{\beta}$$에 대해 OLS 추정량 $$\hat{\boldsymbol{\beta}}_{OLS}$$와 IV 추정량 $$\hat{\boldsymbol{\beta}}_{IV}$$를 도출한다.

**②** 두 추정량 벡터의 차이를 공분산 행렬로 표준화하여, 다음의 하우스만 검정 통계량 $$H$$를 도출한다.

$$
H = (\hat{\boldsymbol{\beta}}_{IV} - \hat{\boldsymbol{\beta}}_{OLS})^\mathsf{T} [\text{Var}(\hat{\boldsymbol{\beta}}_{IV}) - \text{Var}(\hat{\boldsymbol{\beta}}_{OLS})]^{-1} (\hat{\boldsymbol{\beta}}_{IV} - \hat{\boldsymbol{\beta}}_{OLS})
$$

**③** 귀무가설이 참일 때, 이 $$H$$ 통계량은 점근적으로 내생성이 의심되는 설명변수의 개수 $$m$$을 자유도로 하는 카이제곱 분포를 따른다.

$$
H \sim_{asy} \chi^2_m
$$

**④** 계산된 $$H$$ 값이 카이제곱 임계치를 초과하면 귀무가설을 기각하여 내생성이 존재한다고 판단한다.

### 사간-한센 검정 (Sargan-Hansen Test)

도구변수의 수가 내생변수의 수보다 많은 과대식별(overidentified) 상황에서, 사용된 도구변수들이 오차항과 상관관계가 없는 타당한 외생변수인지 검정한다.

- **귀무가설($$H_0$$):** $$\mathbb{E}[\mathbf{Z}^\mathsf{T} u] = \mathbf{0}$$ (모든 도구변수가 오차항과 무상관이다. 즉, 도구변수들이 모두 타당한 외생성을 가진다.)
- **대립가설($$H_1$$):** $$\mathbb{E}[\mathbf{Z}^\mathsf{T} u] \neq \mathbf{0}$$ (적어도 하나 이상의 도구변수가 오차항과 상관되어 있어 타당하지 않다.)

**검정 절차:**
내생변수가 $$k$$개이고 도구변수가 $$m$$개 ($$m > k$$)라고 가정한다.

**①** 도구변수를 이용해 2SLS(또는 IV) 추정을 수행하고 잔차 $$\hat{u}^{IV}_i$$를 추출한다.

**②** 이 잔차를 종속변수로, 모형의 **모든 도구변수와 외생변수**를 설명변수로 하여 보조 회귀식을 추정한다.

$$
\hat{u}^{IV}_i = \gamma_1 Z_{1i} + \gamma_2 Z_{2i} + \cdots + \gamma_m Z_{mi} + \nu_i
$$

**③** 이 보조 회귀식으로부터 결정계수 $$R^2$$를 구한다.

**④** 귀무가설 하에서 다음의 검정 통계량은 점근적으로 자유도가 $$(m-k)$$인 카이제곱 분포를 따른다.

$$
n R^2 \sim_{asy} \chi^2_{m-k}
$$

**⑤** 계산된 통계량이 임계치를 초과하면 귀무가설을 기각하여, 도구변수 세트 중 일부가 타당하지 않다고 결론 내린다.

## Remedy for OVP

OVP를 해결하기 위한 방법으로는 다음이 있으며, 이들에 대해서는 앞으로 살펴볼 것이다.

**①** Control Variable 추가

**②** Proxy Variable 사용

**③** IV 사용

**④** Fixed Effect Model 사용

## 종속변수 관측오차

진정한 데이터 생성 과정(True DGP)이 다음과 같다고 가정한다.

$$
y_i^* = \beta_0 + \beta_1 x_i + u_i
\tag{1}
$$

여기서 고전적 선형회귀 가정인 $$\mathbb{E}(u_i|x_i)=0$$, $$Var(u_i|x_i)=\sigma_u^2$$, $$Cov(u_i, u_j)=0 \ (i \neq j)$$이 성립한다.
우리가 실제로 관측하는 종속변수 $$y_i$$는 진정한 값 $$y_i^*$$에 측정오차 $$e_i$$가 존재한다.

$$
y_i = y_i^* + e_i
\tag{2}
$$

관측된 $$y_i$$를 식 (1)과 (2)를 이용해 원래 회귀식에 대입하여 추정 모형을 구성하면 다음과 같다.

$$
y_i = \beta_0 + \beta_1 x_i + (u_i + e_i) = \beta_0 + \beta_1 x_i + v_i \quad \text{where} \quad v_i \equiv u_i + e_i
\tag{3}
$$

#### 불편성

OLS 기울기 추정량 $$\hat{\beta}_1$$은 다음과 같다.

$$
\hat{\beta}_1 = \frac{\sum_{i=1}^n (x_i - \bar{x})(y_i - \bar{y})}{\sum_{i=1}^n (x_i - \bar{x})^2} = \beta_1 + \frac{\sum_{i=1}^n (x_i - \bar{x})v_i}{\sum_{i=1}^n (x_i - \bar{x})^2}
$$

설명변수 $$X$$에 대한 조건부 기댓값을 취하면,

$$
\mathbb{E}(\hat{\beta}_1 | X) = \beta_1 + \frac{\sum_{i=1}^n (x_i - \bar{x}) \mathbb{E}(v_i | X)}{\sum_{i=1}^n (x_i - \bar{x})^2}
$$

복합오차항의 조건부 기댓값은 가정에 의해 분해되어 소거된다.

$$
\mathbb{E}(v_i | X) = \mathbb{E}(u_i + e_i | X) = \mathbb{E}(u_i | X) + \mathbb{E}(e_i | X) = 0 + 0 = 0 \Longrightarrow \mathbb{E}(\hat{\beta}_1 | X) = \beta_1
$$

고전적 가정을 만족하는 종속변수 측정오차는 **OLS 추정량의 불편성을 훼손하지 않는다**.

#### 일치성

확률수렴과 Slutsky 정리를 적용한다.

$$
\text{plim } \hat{\beta}_1 = \beta_1 + \frac{\text{plim } \frac{1}{n}\sum_{i=1}^n (x_i - \bar{x})v_i}{\text{plim } \frac{1}{n}\sum_{i=1}^n (x_i - \bar{x})^2} = \beta_1 + \frac{\text{Cov}(x, v)}{\text{Var}(x)}
$$

분자의 공분산은 다음과 같이 전개된다.

$$
\text{Cov}(x, v) = \text{Cov}(x, u + e) = \text{Cov}(x, u) + \text{Cov}(x, e) = 0 + 0 = 0
$$

표본 크기 $$n \to \infty$$일 때 $$\beta_1$$으로 수렴하므로 **일치추정량의 성질을 유지한다.**

#### 분산 및 효율성

&#8203;$$u_i$$와 $$e_i$$가 상호 독립이므로 복합오차항의 조건부 분산은 다음과 같다.

$$
Var(v_i | X) = Var(u_i + e_i | X) = Var(u_i|X) + Var(e_i|X) = \sigma_u^2 + \sigma_e^2 \equiv \sigma_v^2
$$

따라서 $$\hat{\beta}_1$$의 분산은 다음과 같이 도출된다.

$$
Var(\hat{\beta}_1 | X) = \frac{\sigma_v^2}{\sum_{i=1}^n (x_i - \bar{x})^2} = \frac{\sigma_u^2 + \sigma_e^2}{\sum_{i=1}^n (x_i - \bar{x})^2}
\tag{4}
$$

측정오차가 없을 때의 분산 $$\left[\frac{\sigma_u^2}{\sum_{i=1}^n (x_i - \bar{x})^2}\right]$$과 비교하면 항상 다음이 성립한다.

$$
Var(\hat{\beta}_1 | X)_{\text{with error}} > Var(\hat{\beta}_1 | X)_{\text{no error}}
\tag{5}
$$

추정량의 중심 위치(기댓값)는 변하지 않으나, **분산이 $$\sigma_e^2 / \sum(x_i-\bar{x**)^2$$ 만큼 팽창하여 OLS의 효율성이 상실된다.}

### 관련 이슈

#### 독립변수 측정오차(CEV)와의 비대칭성

- **종속변수 오차($$y$$):** 설명했듯 불편·일치성이 유지되고 **효율성만 저하**된다.
- **독립변수 오차($$x$$):** 관측값 $$x = x^* + e_x$$ 일 때, OLS 추정량은 $$\text{plim }\hat{\beta}_1 = \beta_1 \left(\frac{\sigma_{x^*}^2}{\sigma_{x^*}^2 + \sigma_{e_x}^2}\right)$$가 되어 0 방향으로 왜곡되는 **감쇠 편의**와 불일치성을 유발한다.

#### 비고전적 측정오차로 인한 내생성

만약 측정오차 $$e_i$$가 설명변수 $$x_i$$와 체계적 상관관계를 가진다면($$Cov(x, e) \neq 0$$), 직교성 가정이 깨진다.

$$
\text{Cov}(x, v) = \text{Cov}(x, u) + \text{Cov}(x, e) = 0 + \text{Cov}(x, e) \neq 0
\tag{6}
$$

- **실증 예시:** 가구 소득($$x$$)이 높을수록 정부 보조금 수급액($$y$$)을 축소 보고하는 경향이 있는 데이터.
- **결과:** $$\mathbb{E}(\hat{\beta}_1) \neq \beta_1$$이 되어 **OLS는 편의를 가지며 불일치**하게 된다. 이 경우 종속변수 오차임에도 도구변수(IV) 추정법이나 패널 GMM 보정이 필수적이다.

#### 측정오차의 평균이 0이 아닌 경우 (절편항 편의)

설문조사의 체계적 과장/축소 등으로 오차의 평균이 0이 아닐 때 모형은 다음과 같이 재편된다.

$$
y_i = (\beta_0 + \mu_e) + \beta_1 x_i + (u_i + e_i - \mu_e)
\tag{7}
$$

- 기울기 $$\hat{\beta}_1$$은 여전히 불편·일치성을 유지한다.
- 그러나 절편항 추정량은 $$\text{plim }\hat{\beta}_0 = \beta_0 + \mu_e$$가 되어 **절편항에 상향 혹은 하향 편의가 발생**한다. 정책의 한계효과(slope) 분석에는 문제가 없으나, 수준(level)을 평가하거나 예측값 $$\hat{y}$$을 도출할 때는 왜곡이 발생한다.

#### 통계적 검정력(Statistical Power) 약화와 제2종 오류

분산 팽창은 추정된 표준오차 $$\text{se}(\hat{\beta}_1)$$의 과대평가로 이어진다.

- $$t$$-통계량 $$\left(\frac{\hat{\beta}_1}{\text{se}(\hat{\beta}_1)}\right)$$의 분모가 커지므로 통계량의 절대값이 작아진다.
- 이는 귀무가설($$H_0: \beta_1 = 0$$)이 거짓임에도 불구하고 이를 기각하지 못하는 **제2종 오류(Type II Error)의 발생 확률을 높인다.**
- **연구자의 대응:** 노이즈가 심한 종속변수를 다룰 때는 표본 크기 $$n$$을 대폭 늘려 분산 식의 분모항 $$\sum(x_i-\bar{x})^2$$을 키움으로써 정밀도 손실을 보상해야 한다.

#### 가법(Additive) 오차 vs 승법(Multiplicative) 오차

- 관측값이 $$y_i = y_i^* \cdot \eta_i$$ 형태의 승법 오차일 경우, 선형 모형에 그대로 OLS를 적용하면 불편성이 깨진다.
- 단, 원래 모형이 로그-선형(Log-linear: $$\ln y^* = \beta_0 + \beta_1 x + u$$) 형태라면 양변 로그 취함에 따라 $$\ln y = \ln y^* + \ln \eta$$가 되어 표준적인 가법 측정오차 모형으로 자연 치환된다. 따라서 **종속변수 측정오차가 의심될 때 로그 변환은 모형의 강건성을 높이는 좋은 수단이 된다.**
