---
title: "[Econometrics] Partitioned Regression"
date: 2026-11-18 16:20:00 +0900
categories: [Econometrics]
tags: [OLS, FWL, 분할회귀, 잔차생성행렬]
math: true
toc: true
---

## Introduction

프리쉬-워-로벨 정리(FWL Thm)는 다중회귀분석에서 특정 독립변수의 순수한 효과를 다른 변수들의 영향을 통제한 상태로 도출하는 과정을 보여준다.

## Scalar Version

독립변수가 2개인 단순한 다중회귀모형을 가정한다.

$$
Y_i = \beta_0 + \beta_1 X_{1i} + \beta_2 X_{2i} + u_i
$$

목표는 $$X_{2i}$$의 영향을 통제한 상태에서 $$X_{1i}$$가 $$Y_i$$에 미치는 순수한 효과인 $$\hat{\beta}_1$$을 구하는 것이다.

**①** $$X_{1i}$$를 종속변수로, $$X_{2i}$$를 독립변수로 두고 단순회귀분석을 실행하여 잔차 $$\hat{r}_{1i}$$를 구한다.

$$
X_{1i} = \hat{\gamma}_0 + \hat{\gamma}_1 X_{2i} + \hat{r}_{1i}
$$

여기서 $$\hat{r}_{1i}$$는 $$X_{2i}$$로는 설명되지 않는 $$X_{1i}$$만의 고유한 정보(순수한 변동)이다.

**②** 원래의 종속변수 $$Y_i$$를 방금 구한 잔차 $$\hat{r}_{1i}$$에 대해 단순회귀분석을 실행한다.

**③** 이때 도출되는 기울기 계수가 다중회귀분석에서의 $$\hat{\beta}_1$$과 정확히 일치한다.

$$
\hat{\beta}_1 = \frac{\sum_{i=1}^n \hat{r}_{1i} Y_i}{\sum_{i=1}^n \hat{r}_{1i}^2}
$$

독립변수가 3개인 다중회귀모형으로 확장해보자.

$$
Y_i = \beta_0 + \beta_1 X_{1i} + \beta_2 X_{2i} + \beta_3 X_{3i} + u_i
$$

목표는 $$X_{2i}$$와 $$X_{3i}$$의 영향을 모두 통제한 상태에서 $$X_{1i}$$가 $$Y_i$$에 미치는 순수한 효과인 $$\hat{\beta}_1$$을 구하는 것이다.

**①** $$X_{1i}$$를 종속변수로, 통제하고자 하는 나머지 변수인 $$X_{2i}$$와 $$X_{3i}$$를 독립변수로 두고 다중회귀분석을 실행하여 잔차 $$\hat{r}_{1i}$$를 구한다.

$$
X_{1i} = \hat{\alpha}_0 + \hat{\alpha}_1 X_{2i} + \hat{\alpha}_2 X_{3i} + \hat{r}_{1i}
$$

여기서 도출된 잔차 $$\hat{r}_{1i}$$는 $$X_{2i}$$와 $$X_{3i}$$로는 전혀 설명되지 않는 $$X_{1i}$$만의 순수한 고유 변동을 의미한다.

**②** 원래의 종속변수 $$Y_i$$를 방금 구한 잔차 $$\hat{r}_{1i}$$에 대해 단순회귀분석을 실행한다.<br>(참고로 종속변수 $$Y_i$$ 역시 $$X_{2i}$$와 $$X_{3i}$$에 대해 회귀하여 얻은 잔차 $$\hat{r}_{Yi}$$를 도출한 후, 이를 $$\hat{r}_{1i}$$에 대해 회귀분석하여도 완벽히 동일한 결과를 얻는다.)

**③** 이때 도출되는 기울기 계수가 원 모형의 $$\hat{\beta}_1$$과 정확히 일치한다.

$$
\hat{\beta}_1 = \frac{\sum_{i=1}^n \hat{r}_{1i} Y_i}{\sum_{i=1}^n \hat{r}_{1i}^2}
$$

## Matrix Version

다수의 변수가 포함된 일반적인 다중회귀모형 $$Y = X\beta + u$$를 두 그룹의 변수 행렬 $$X_1$$과 $$X_2$$로 Partition한다.

$$
Y = X_1\beta_1 + X_2\beta_2 + u
$$

OLS의 Normal Equations $$X^\mathsf{T}X\hat{\beta} = X^\mathsf{T}Y$$를 분할된 행렬 형태로 쓰면 다음과 같다.

$$
\begin{bmatrix} X_1^\mathsf{T}X_1 & X_1^\mathsf{T}X_2 \\ X_2^\mathsf{T}X_1 & X_2^\mathsf{T}X_2 \end{bmatrix} \begin{bmatrix} \hat{\beta}_1 \\ \hat{\beta}_2 \end{bmatrix} = \begin{bmatrix} X_1^\mathsf{T}Y \\ X_2^\mathsf{T}Y \end{bmatrix}
$$

### 도출 과정

**①** 두 번째 줄의 방정식을 $$\hat{\beta}_2$$에 대해 정리한다.

$$
X_2^\mathsf{T}X_1 \hat{\beta}_1 + X_2^\mathsf{T}X_2 \hat{\beta}_2 = X_2^\mathsf{T}Y
$$

$$
\hat{\beta}_2 = (X_2^\mathsf{T}X_2)^{-1}(X_2^\mathsf{T}Y - X_2^\mathsf{T}X_1 \hat{\beta}_1)
$$

**②** 구한 $$\hat{\beta}_2$$를 첫 번째 줄의 방정식에 대입한다.

$$
X_1^\mathsf{T}X_1 \hat{\beta}_1 + X_1^\mathsf{T}X_2 [ (X_2^\mathsf{T}X_2)^{-1}(X_2^\mathsf{T}Y - X_2^\mathsf{T}X_1 \hat{\beta}_1) ] = X_1^\mathsf{T}Y
$$

**③** 식을 전개하고 $$\hat{\beta}_1$$로 묶는다.

$$
X_1^\mathsf{T}X_1 \hat{\beta}_1 - X_1^\mathsf{T}X_2(X_2^\mathsf{T}X_2)^{-1}X_2^\mathsf{T}X_1 \hat{\beta}_1 = X_1^\mathsf{T}Y - X_1^\mathsf{T}X_2(X_2^\mathsf{T}X_2)^{-1}X_2^\mathsf{T}Y
$$

$$
X_1^\mathsf{T}[I - X_2(X_2^\mathsf{T}X_2)^{-1}X_2^\mathsf{T}]X_1 \hat{\beta}_1 = X_1^\mathsf{T}[I - X_2(X_2^\mathsf{T}X_2)^{-1}X_2^\mathsf{T}]Y
$$

**④** 잔차 생성 행렬(Residual Maker Matrix) $$M_2$$의 도입.

$$
M_2 = I - X_2(X_2^\mathsf{T}X_2)^{-1}X_2^\mathsf{T}
$$

이 행렬은 $$X_2$$로 설명되는 부분을 제거하고 잔차만 남기는 projection matrix이다. 이를 식에 대입하면 다음과 같다.

$$
X_1^\mathsf{T}M_2 X_1 \hat{\beta}_1 = X_1^\mathsf{T}M_2 Y
$$

**⑤** 최종적으로 양변에 역행렬을 곱하여 $$\hat{\beta}_1$$을 도출한다.

$$
\hat{\beta}_1 = (X_1^\mathsf{T}M_2 X_1)^{-1} X_1^\mathsf{T}M_2 Y
$$

### 결론

행렬 $$M_2$$는 대칭행렬($$M_2 = M_2^\mathsf{T}$$)이자 idempotent matrix ($$M_2 M_2 = M_2$$)이다. 이를 활용해 최종 식을 다시 쓰면 다음과 같다.

$$
\hat{\beta}_1 = (X_1^\mathsf{T} M_2^\mathsf{T} M_2 X_1)^{-1} (X_1^\mathsf{T} M_2^\mathsf{T} M_2 Y) = ({X_1^*}^\mathsf{T} X_1^*)^{-1} {X_1^*}^\mathsf{T} Y^*
$$

단, $$X_1^* = M_2 X_1$$, $$Y^* = M_2 Y$$이다. 이는 $$Y$$와 $$X_1$$에서 $$X_2$$의 영향을 모두 제거해버린 잔차 행렬끼리 회귀분석을 돌린 결과와 정확히 동일함을 보여준다.
