
# WinForm статистический анализ и проверка гипотез

[![C++](https://img.shields.io/badge/C++-CLI/Windows-blue.svg)](https://docs.microsoft.com/ru-ru/cpp/dotnet/dotnet-programming-with-cpp-cli-visual-cpp)
[![Статус](https://img.shields.io/badge/Статус-Завершено-brightgreen.svg)]()
![GitHub Repo stars](https://img.shields.io/github/stars/nedashh/Statistical-analysis)
![GitHub forks](https://img.shields.io/github/forks/nedashh/Statistical-analysis)

**Русский**

Десктопное приложение для Windows (C++/CLI), предназначенное для статистического оценивания и проверки гипотез. Разработано в рамках курса **«Алгоритмы и обработка данных»**.

---

## 📋 Обзор

Проект реализует широкий набор методов математической статистики для анализа экспериментальных данных, включая оценивание параметров распределений (ММП, МНК) и проверку статистических гипотез (параметрические и непараметрические критерии).

### Статистическая оценка
- **Метод максимального правдоподобия (ММП)**
  - Нормальное распределение
  - Распределение Вейбулла
- **Метод наименьших квадратов (МНК)**
  - Нормальное распределение
  - Распределение Вейбулла

### Параметрическая проверка гипотез
- **Критерий аномальности**
  - Критерий Граббса
- **Критерий однородности**
  - Критерий Стьюдента
  - Однофакторный дисперсионный анализ
  - Критерий Бартлетта
- **Критерий согласия**
  - Критерий Шапиро–Уилка

### Непараметрические проверки гипотез
- Критерий Краскела-Уоллиса
- Критерий знаковых рангов Уилкоксона
- Двухвыборочный критерий Уилкоксона


---

## 🚀 Начало работы

### Требования
- Visual Studio 2017 или новее
- .NET Framework 4.5+
- Библиотеки Boost C++

### Установка

```bash
git clone https://github.com/nedashh/Statistical-analysis.git
cd Statistical-analysis
```

1. Установите Boost C++ Libraries  
2. Установите Visual Studio с поддержкой C++/CLI  
3. Откройте `Methods.sln`  
4. Укажите пути к библиотекам Boost  
5. Выберите платформу (x86 / x64)  
6. Соберите проект (рекомендуется Release)

---

# WinForm Statistical Analysis and Hypothesis Testing

**English** 

A Windows desktop application (C++/CLI) designed for statistical parameter estimation and hypothesis testing. Developed as part of the **“Algorithms and Data Processing”** course.

---

## 📋 Overview

The project implements a wide range of mathematical statistics methods for the analysis of experimental data, including distribution parameter estimation (Maximum Likelihood Estimation and Least Squares Method) and statistical hypothesis testing (parametric and non-parametric tests).

### Statistical Estimation
- **Maximum Likelihood Estimation (MLE)**
  - Normal distribution
  - Weibull distribution
- **Least Squares Method (LSM)**
  - Normal distribution
  - Weibull distribution

### Parametric Hypothesis Tests
- **Anomaly Detection**
  - Grubbs test for outliers
- **Homogeneity Tests**
  - Student’s t-test
  - One-way ANOVA
  - Bartlett’s test
- **Goodness-of-Fit**
  - Shapiro–Wilk normality test

### Non-Parametric Hypothesis Tests
- Wilcoxon signed-rank test
- Mann–Whitney U test
- Kruskal–Wallis test

---

## 🚀 Getting Started

### Requirements
- Visual Studio 2017 or later
- .NET Framework 4.5+
- Boost C++ Libraries

### Installation

```bash
git clone https://github.com/nedashh/Statistical-analysis.git
cd Statistical-analysis
Install Boost C++ Libraries
```

1. Install Visual Studio with C++/CLI support
2. Open `Methods.sln` 
3. Configure Boost include and library paths
4. Select the target platform (x86 / x64)
5. Build the project (Release configuration recommended)
