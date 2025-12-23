#pragma once
#include "DataGenerator.h"
#include "stat.h"
#include <msclr/marshal_cppstd.h>
#include <string>

using namespace System;
using namespace System::ComponentModel;
using namespace System::Collections;
using namespace System::Collections::Generic;
using namespace System::Windows::Forms;
using namespace System::Data;
using namespace System::Drawing;
using namespace System::IO;
using namespace System::Diagnostics;
using namespace System::Windows::Forms::DataVisualization::Charting;
using namespace System::Globalization;

public ref class MainForm : public System::Windows::Forms::Form
{
public:
    MainForm(void)
    {
        InitializeComponent();
        InitializeTreeView();
    }

private:
    System::Windows::Forms::TreeView^ treeView;
    System::Windows::Forms::SplitContainer^ splitContainer;
    System::Windows::Forms::ToolStrip^ toolStrip;
    System::Windows::Forms::StatusStrip^ statusStrip;
    System::Windows::Forms::RichTextBox^ outputBox;
    Chart^ chart;
    System::ComponentModel::Container^ components;

    ToolStripButton^ btnLoadInp;
    ToolStripButton^ btnLoadOut;
    ToolStripButton^ btnRunTest;
    ToolStripButton^ btnClear;
    ToolStripButton^ btnPlot;

    String^ currentMethod;

#pragma region Windows Form Designer generated code

    void InitializeComponent(void)
    {
        this->components = (gcnew System::ComponentModel::Container());
        this->treeView = (gcnew System::Windows::Forms::TreeView());
        this->splitContainer = (gcnew System::Windows::Forms::SplitContainer());
        this->outputBox = (gcnew System::Windows::Forms::RichTextBox());
        this->chart = (gcnew System::Windows::Forms::DataVisualization::Charting::Chart());
        this->toolStrip = (gcnew System::Windows::Forms::ToolStrip());
        this->btnLoadInp = (gcnew System::Windows::Forms::ToolStripButton());
        this->btnLoadOut = (gcnew System::Windows::Forms::ToolStripButton());
        this->btnRunTest = (gcnew System::Windows::Forms::ToolStripButton());
        this->btnClear = (gcnew System::Windows::Forms::ToolStripButton());
        this->btnPlot = (gcnew System::Windows::Forms::ToolStripButton());
        this->statusStrip = (gcnew System::Windows::Forms::StatusStrip());

        this->splitContainer->SuspendLayout();
        this->toolStrip->SuspendLayout();
        this->SuspendLayout();

        this->treeView->Dock = System::Windows::Forms::DockStyle::Left;
        this->treeView->Location = System::Drawing::Point(0, 25);
        this->treeView->Name = L"treeView";
        this->treeView->Size = System::Drawing::Size(250, 616);
        this->treeView->TabIndex = 0;
        this->treeView->AfterSelect += gcnew System::Windows::Forms::TreeViewEventHandler(this, &MainForm::treeView_AfterSelect);

        this->splitContainer->Dock = System::Windows::Forms::DockStyle::Fill;
        this->splitContainer->Location = System::Drawing::Point(250, 25);
        this->splitContainer->Name = L"splitContainer";
        this->splitContainer->Orientation = System::Windows::Forms::Orientation::Horizontal;
        this->splitContainer->Size = System::Drawing::Size(734, 616);
        this->splitContainer->SplitterDistance = 300;
        this->splitContainer->TabIndex = 1;

        this->outputBox->Dock = System::Windows::Forms::DockStyle::Fill;
        this->outputBox->Location = System::Drawing::Point(0, 0);
        this->outputBox->Name = L"outputBox";
        this->outputBox->Size = System::Drawing::Size(734, 300);
        this->outputBox->TabIndex = 0;

        this->chart->Dock = System::Windows::Forms::DockStyle::Fill;
        this->chart->Location = System::Drawing::Point(0, 0);
        this->chart->Name = L"chart";
        this->chart->Size = System::Drawing::Size(734, 312);
        this->chart->TabIndex = 1;
        this->chart->Text = L"График";

        ChartArea^ chartArea = gcnew ChartArea();
        chartArea->Name = "MainArea";
        chartArea->AxisX->Title = "X";
        chartArea->AxisY->Title = "Y";
        this->chart->ChartAreas->Add(chartArea);

        Legend^ legend = gcnew Legend();
        this->chart->Legends->Add(legend);

        this->splitContainer->Panel1->Controls->Add(this->outputBox);
        this->splitContainer->Panel2->Controls->Add(this->chart);

        this->toolStrip->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(5) {
            this->btnLoadInp, this->btnLoadOut, this->btnRunTest, this->btnClear, this->btnPlot
        });
        this->toolStrip->Location = System::Drawing::Point(0, 0);
        this->toolStrip->Name = L"toolStrip";
        this->toolStrip->Size = System::Drawing::Size(984, 25);
        this->toolStrip->TabIndex = 2;

        this->btnLoadInp->Text = L"Загрузить INP";
        this->btnLoadInp->Click += gcnew System::EventHandler(this, &MainForm::btnLoadInp_Click);

        this->btnLoadOut->Text = L"Загрузить OUT";
        this->btnLoadOut->Click += gcnew System::EventHandler(this, &MainForm::btnLoadOut_Click);

        this->btnRunTest->Text = L"Выполнить тест";
        this->btnRunTest->Click += gcnew System::EventHandler(this, &MainForm::btnRunTest_Click);

        this->btnClear->Text = L"Очистить";
        this->btnClear->Click += gcnew System::EventHandler(this, &MainForm::btnClear_Click);

        this->btnPlot->Text = L"Построить график";
        this->btnPlot->Click += gcnew System::EventHandler(this, &MainForm::btnPlot_Click);

        this->statusStrip->Location = System::Drawing::Point(0, 641);
        this->statusStrip->Name = L"statusStrip";
        this->statusStrip->Size = System::Drawing::Size(984, 22);
        this->statusStrip->TabIndex = 3;

        this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
        this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
        this->ClientSize = System::Drawing::Size(984, 663);
        this->Controls->Add(this->splitContainer);
        this->Controls->Add(this->treeView);
        this->Controls->Add(this->toolStrip);
        this->Controls->Add(this->statusStrip);
        this->Name = L"MainForm";
        this->Text = L"Статистический анализ";
        this->StartPosition = System::Windows::Forms::FormStartPosition::CenterScreen;

        this->splitContainer->ResumeLayout(false);
        this->toolStrip->ResumeLayout(false);
        this->toolStrip->PerformLayout();
        this->ResumeLayout(false);
        this->PerformLayout();
    }

#pragma endregion

private:
    void InitializeTreeView()
    {
        treeView->Nodes->Clear();

        TreeNode^ root = gcnew TreeNode("Статистический анализ");

        TreeNode^ estimationNode = gcnew TreeNode("Статистическое оценивание");
        TreeNode^ parametricNode = gcnew TreeNode("Параметрические методы");

        TreeNode^ mleNode = gcnew TreeNode("Метод максимального правдоподобия");
        mleNode->Nodes->Add("Нормальный закон")->Tag = "MLE_Normal";
        mleNode->Nodes->Add("Закон Вейбулла")->Tag = "MLE_Weibull";

        TreeNode^ mlsNode = gcnew TreeNode("Метод наименьших квадратов");
        mlsNode->Nodes->Add("Нормальный закон")->Tag = "MLS_Normal";
        mlsNode->Nodes->Add("Закон Вейбулла")->Tag = "MLS_Weibull";

        parametricNode->Nodes->Add(mleNode);
        parametricNode->Nodes->Add(mlsNode);
        estimationNode->Nodes->Add(parametricNode);

        TreeNode^ hypothesesNode = gcnew TreeNode("Проверка статистических гипотез");

        // PARAMETRIC CRITERIA
        TreeNode^ parametricTestsNode = gcnew TreeNode("Параметрические критерии");

        TreeNode^ anomalyNode = gcnew TreeNode("Критерий аномальности");
        anomalyNode->Nodes->Add("Критерий Граббса")->Tag = "GrubbsTest";

        TreeNode^ homogeneityNode = gcnew TreeNode("Критерий однородности");
        homogeneityNode->Nodes->Add("Критерий Стьюдента")->Tag = "StudentTest";
        homogeneityNode->Nodes->Add("Однофакторный дисперсионный анализ")->Tag = "ANOVA";
        homogeneityNode->Nodes->Add("Критерий Бартлетта")->Tag = "BartlettTest";

        // ACCEPTANCE CRITERIA (parametric)
        TreeNode^ goodnessOfFitNode = gcnew TreeNode("Критерии согласия");
        goodnessOfFitNode->Nodes->Add("Критерий Шапиро-Уилка")->Tag = "ShapiroWilkTest";

        parametricTestsNode->Nodes->Add(anomalyNode);
        parametricTestsNode->Nodes->Add(homogeneityNode);
        parametricTestsNode->Nodes->Add(goodnessOfFitNode);

        // NONPARAMETRIC CRITERIA
        TreeNode^ nonParametricTestsNode = gcnew TreeNode("Непараметрические критерии");
        TreeNode^ wilcoxonNode = gcnew TreeNode("Критерии Уилкоксона");
        wilcoxonNode->Nodes->Add("Критерий знаковых рангов Уилкоксона")->Tag = "WilcoxonSignedRankTest";
        wilcoxonNode->Nodes->Add("Двухвыборочный критерий Уилкоксона")->Tag = "WilcoxonRankSumTest";
        nonParametricTestsNode->Nodes->Add("Критерий Краскела-Уоллиса")->Tag = "KruskalWallisTest";
        nonParametricTestsNode->Nodes->Add(wilcoxonNode);

        hypothesesNode->Nodes->Add(parametricTestsNode);
        hypothesesNode->Nodes->Add(nonParametricTestsNode);

        root->Nodes->Add(estimationNode);
        root->Nodes->Add(hypothesesNode);


        treeView->Nodes->Add(root);
        treeView->ExpandAll();
    }

    void treeView_AfterSelect(Object^ sender, TreeViewEventArgs^ e)
    {
        if (e->Node->Tag != nullptr)
        {
            currentMethod = e->Node->Tag->ToString();
            UpdateStatus("Выбран метод: " + currentMethod);
        }
    }

    void btnLoadInp_Click(Object^ sender, EventArgs^ e)
    {
        if (String::IsNullOrEmpty(currentMethod))
        {
            MessageBox::Show("Сначала выберите метод из дерева", "Ошибка");
            return;
        }

        try
        {
            Directory::CreateDirectory("Inp");
            String^ inpFile = "Inp/" + currentMethod + ".inp";
            std::string content;

            if (currentMethod == "MLE_Normal") content = DataGenerator::generateMLE_Normal();
            else if (currentMethod == "MLE_Weibull") content = DataGenerator::generateMLE_Weibull();
            else if (currentMethod == "MLS_Normal") content = DataGenerator::generateMLS_Normal();
            else if (currentMethod == "MLS_Weibull") content = DataGenerator::generateMLS_Weibull();
            else if (currentMethod == "GrubbsTest") content = DataGenerator::generateGrubbsTest();
            else if (currentMethod == "StudentTest") content = DataGenerator::generateStudentTest();
            else if (currentMethod == "BartlettTest") content = DataGenerator::generateBartlettTest();
            else if (currentMethod == "ANOVA") content = DataGenerator::generateANOVA();
            else if (currentMethod == "ShapiroWilkTest") content = DataGenerator::generateShapiroWilkTest();
            else if (currentMethod == "WilcoxonSignedRankTest") content = DataGenerator::generateWilcoxonSignedRankTest();
            else if (currentMethod == "WilcoxonRankSumTest") content = DataGenerator::generateWilcoxonRankSumTest();
            else if (currentMethod == "KruskalWallisTest") content = DataGenerator::generateKruskalWallisTest();
            else
            {
                content = "sample_size\n20\ndata\n";
                for (int i = 0; i < 20; i++) {
                    double value = 10.0 + DataGenerator::randomDouble(0, 5.0);
                    content += std::to_string(value).substr(0, 5) + " ";
                }
            }

            String^ managedContent = gcnew String(content.c_str());
            File::WriteAllText(inpFile, managedContent);
            outputBox->Text = "=== ВХОДНЫЕ ДАННЫЕ ===\nМетод: " + currentMethod + "\n\n" + managedContent;
            UpdateStatus("Создан файл: " + inpFile);
        }
        catch (Exception^ ex)
        {
            MessageBox::Show("Ошибка: " + ex->Message, "Ошибка");
        }
    }

    void btnLoadOut_Click(Object^ sender, EventArgs^ e)
    {
        if (String::IsNullOrEmpty(currentMethod))
        {
            MessageBox::Show("Сначала выберите метод из дерева", "Ошибка");
            return;
        }

        try
        {
            String^ outFile = "Out/" + currentMethod + ".out";
            if (File::Exists(outFile))
            {
                String^ content = File::ReadAllText(outFile);
                outputBox->Text = "=== РЕЗУЛЬТАТЫ ===\nМетод: " + currentMethod + "\n\n" + content;
                UpdateStatus("Загружены результаты: " + currentMethod);
            }
            else
            {
                MessageBox::Show("Файл результатов не найден. Сначала выполните тест.", "Информация");
            }
        }
        catch (Exception^ ex)
        {
            MessageBox::Show("Ошибка: " + ex->Message, "Ошибка");
        }
    }

    void btnRunTest_Click(Object^ sender, EventArgs^ e)
    {
        if (String::IsNullOrEmpty(currentMethod))
        {
            MessageBox::Show("Сначала выберите метод из дерева", "Ошибка");
            return;
        }

        try
        {
            String^ inpFile = "Inp/" + currentMethod + ".inp";
            if (!File::Exists(inpFile))
            {
                MessageBox::Show("Сначала создайте входные данные через 'Загрузить INP'", "Ошибка");
                return;
            }

            Directory::CreateDirectory("Out");

            outputBox->Text = "Выполняется " + currentMethod + "...";
            Application::DoEvents();

            bool testExecuted = false;

            if (currentMethod == "MLE_Normal") {
                MLE_Normal("MLE_Normal");
                testExecuted = true;
            }
            else if (currentMethod == "MLE_Weibull") {
                MLE_Weibull("MLE_Weibull");
                testExecuted = true;
            }
            else if (currentMethod == "MLS_Normal") {
                MLS_Normal("MLS_Normal");
                testExecuted = true;
            }
            else if (currentMethod == "MLS_Weibull") {
                MLS_Weibull("MLS_Weibull");
                testExecuted = true;
            }
            else if (currentMethod == "GrubbsTest") {
                GrubbsTest("GrubbsTest");
                testExecuted = true;
            }
            else if (currentMethod == "StudentTest") {
                StudentTest("StudentTest");
                testExecuted = true;
            }
            else if (currentMethod == "BartlettTest") {
                BartlettTest("BartlettTest");
                testExecuted = true;
            }
            else if (currentMethod == "ANOVA") {
                ANOVA("ANOVA");
                testExecuted = true;
            }
            else if (currentMethod == "ShapiroWilkTest") {
                ShapiroWilkTest("ShapiroWilkTest");
                testExecuted = true;
            }
            else if (currentMethod == "WilcoxonSignedRankTest") {
                WilcoxonSignedRankTest("WilcoxonSignedRankTest");
                testExecuted = true;
            }
            else if (currentMethod == "WilcoxonRankSumTest") {
                WilcoxonRankSumTest("WilcoxonRankSumTest");
                testExecuted = true;
            }
            else if (currentMethod == "KruskalWallisTest") {
                KruskalWallisTest("KruskalWallisTest");
                testExecuted = true;
            }
            if (!testExecuted) {
                MessageBox::Show("Данный метод еще не реализован", "Информация");
                return;
            }

            String^ outFile = "Out/" + currentMethod + ".out";
            if (File::Exists(outFile)) {
                String^ content = File::ReadAllText(outFile);
                outputBox->Text = "=== РЕЗУЛЬТАТЫ ===\nМетод: " + currentMethod + "\n\n" + content;
                UpdateStatus("Тест выполнен: " + currentMethod);
            }
            else {
                outputBox->Text = "ОШИБКА: Файл результатов не создан: " + outFile +
                    "\nПроверьте:\n1. Правильность формата входных данных\n2. Реализацию функции " + currentMethod;
                UpdateStatus("Ошибка выполнения: " + currentMethod);
            }
        }
        catch (Exception^ ex) {
            MessageBox::Show("Ошибка выполнения теста: " + ex->Message, "Ошибка");
        }
    }

    void btnClear_Click(Object^ sender, EventArgs^ e)
    {
        outputBox->Clear();
        chart->Series->Clear();
        UpdateStatus("Очищено");
    }

    void btnPlot_Click(Object^ sender, EventArgs^ e)
    {
        if (String::IsNullOrEmpty(currentMethod))
        {
            MessageBox::Show("Сначала выберите метод из дерева");
            return;
        }

        chart->Series->Clear();
        chart->Legends->Clear();
        chart->Titles->Clear();

        Legend^ legend = gcnew Legend();
        legend->Docking = Docking::Top;
        legend->Alignment = Drawing::StringAlignment::Center;
        chart->Legends->Add(legend);

        chart->Titles->Add("Функция распределения - " + currentMethod);

        try
        {
            if (currentMethod->StartsWith("MLE_") || currentMethod->StartsWith("MLS_"))
            {
                PlotDistributionWithCI(currentMethod);
            }
            else
            {
                PlotSimpleDistribution(currentMethod);
            }
        }
        catch (Exception^ ex)
        {
            MessageBox::Show("Ошибка построения графика: " + ex->Message, "Ошибка");
        }
    }

private:
    void PlotDistributionWithCI(String^ method)
    {
        chart->Series->Clear();

        String^ mainFile = "Out/" + method + "_main.dat";
        String^ lowerFile = "Out/" + method + "_lower.dat";
        String^ upperFile = "Out/" + method + "_upper.dat";

        if (!File::Exists(mainFile))
        {
            MessageBox::Show("Основной файл данных не найден: " + mainFile, "Ошибка");
            return;
        }

        Series^ mainSeries = gcnew Series("Основное распределение");
        mainSeries->ChartType = SeriesChartType::Line;
        mainSeries->Color = Color::Blue;
        mainSeries->BorderWidth = 3;
        mainSeries->MarkerStyle = MarkerStyle::None;

        Series^ lowerSeries = gcnew Series("Нижняя граница");
        lowerSeries->ChartType = SeriesChartType::Line;
        lowerSeries->Color = Color::Red;
        lowerSeries->BorderDashStyle = ChartDashStyle::Dash;
        lowerSeries->BorderWidth = 1;
        lowerSeries->MarkerStyle = MarkerStyle::None;

        Series^ upperSeries = gcnew Series("Верхняя граница");
        upperSeries->ChartType = SeriesChartType::Line;
        upperSeries->Color = Color::Red;
        upperSeries->BorderDashStyle = ChartDashStyle::Dash;
        upperSeries->BorderWidth = 1;
        upperSeries->MarkerStyle = MarkerStyle::None;

        if (ReadDataFromFile(mainFile, mainSeries))
        {
            chart->Series->Add(mainSeries);
            UpdateStatus("Основная серия: " + mainSeries->Points->Count + " точек");
        }

        if (File::Exists(lowerFile) && ReadDataFromFile(lowerFile, lowerSeries))
        {
            chart->Series->Add(lowerSeries);
            UpdateStatus("Нижняя граница: " + lowerSeries->Points->Count + " точек");
        }

        if (File::Exists(upperFile) && ReadDataFromFile(upperFile, upperSeries))
        {
            chart->Series->Add(upperSeries);
            UpdateStatus("Верхняя граница: " + upperSeries->Points->Count + " точек");
        }

        chart->ChartAreas[0]->AxisX->Title = "Вероятность P";
        chart->ChartAreas[0]->AxisY->Title = "Квантиль Xp";
        chart->ChartAreas[0]->AxisX->Minimum = 0;
        chart->ChartAreas[0]->AxisX->Maximum = 1;

        chart->ChartAreas[0]->RecalculateAxesScale();

        UpdateStatus("ГРАФИК ПОСТРОЕН! Всего серий: " + chart->Series->Count);
    }

    bool ReadDataFromFile(String^ filename, Series^ series)
    {
        try
        {
            if (!File::Exists(filename))
                return false;

            cli::array<String^>^ lines = File::ReadAllLines(filename);
            int pointsAdded = 0;

            for each(String ^ line in lines)
            {
                String^ trimmed = line->Trim();
                if (trimmed->StartsWith("#") || String::IsNullOrEmpty(trimmed))
                    continue;

                cli::array<String^>^ parts = trimmed->Split(
                    gcnew cli::array<wchar_t>{' ', '\t'},
                    StringSplitOptions::RemoveEmptyEntries);

                if (parts->Length >= 2)
                {
                    double x, y;
                    if (Double::TryParse(parts[0], NumberStyles::Any, CultureInfo::InvariantCulture, x) &&
                        Double::TryParse(parts[1], NumberStyles::Any, CultureInfo::InvariantCulture, y))
                    {
                        series->Points->AddXY(x, y);
                        pointsAdded++;
                    }
                }
            }

            return (pointsAdded > 0);
        }
        catch (Exception^ ex)
        {
            UpdateStatus("Ошибка чтения " + Path::GetFileName(filename) + ": " + ex->Message);
            return false;
        }
    }

    int GetTotalPoints()
    {
        int total = 0;
        for each(Series ^ series in chart->Series)
        {
            total += series->Points->Count;
        }
        return total;
    }
    void TestChart()
    {
        chart->Series->Clear();

        Series^ testSeries = gcnew Series("Тестовые данные");
        testSeries->ChartType = SeriesChartType::Line;
        testSeries->Color = Color::Red;
        testSeries->BorderWidth = 3;

        testSeries->Points->AddXY(0.1, 1.0);
        testSeries->Points->AddXY(0.5, 2.0);
        testSeries->Points->AddXY(0.9, 3.0);

        chart->Series->Add(testSeries);

        chart->ChartAreas[0]->AxisX->Minimum = 0;
        chart->ChartAreas[0]->AxisX->Maximum = 1;
        chart->ChartAreas[0]->AxisY->Minimum = 0;
        chart->ChartAreas[0]->AxisY->Maximum = 4;

        chart->Invalidate();
        UpdateStatus("Тестовый график построен");
    }

    bool AddDataFromFile(Series^ series, String^ filename)
    {
        try
        {
            if (!File::Exists(filename))
            {
                UpdateStatus("Файл не существует: " + filename);
                return false;
            }

            cli::array<String^>^ allLines = File::ReadAllLines(filename);
            UpdateStatus("Файл " + Path::GetFileName(filename) + ": " + allLines->Length + " строк");

            bool dataFound = false;
            int dataPoints = 0;

            for each(String ^ line in allLines)
            {
                String^ trimmedLine = line->Trim();

                if (trimmedLine->StartsWith("#") || String::IsNullOrEmpty(trimmedLine))
                    continue;

                cli::array<String^>^ parts = trimmedLine->Split(
                    gcnew cli::array<wchar_t>{' ', '\t'},
                    StringSplitOptions::RemoveEmptyEntries);

                if (parts->Length >= 2)
                {
                    double x, y;
                    if (Double::TryParse(parts[0], x) && Double::TryParse(parts[1], y))
                    {
                        series->Points->AddXY(x, y);
                        dataFound = true;
                        dataPoints++;
                    }
                }
            }

            UpdateStatus(Path::GetFileName(filename) + ": " + dataPoints + " точек");

            if (dataPoints > 0)
            {
                String^ pointsInfo = "Первые точки " + Path::GetFileName(filename) + ": ";
                for (int i = 0; i < Math::Min(3, series->Points->Count); i++)
                {
                    pointsInfo += String::Format("({0:F3}, {1:F3}) ",
                        series->Points[i]->XValue, series->Points[i]->YValues[0]);
                }
                UpdateStatus(pointsInfo);
            }

            return dataFound;
        }
        catch (Exception^ ex)
        {
            UpdateStatus("Ошибка: " + Path::GetFileName(filename) + " - " + ex->Message);
            return false;
        }
    }

    void PlotSimpleDistribution(String^ method)
    {
        String^ outFile = "Out/" + method + ".out";
        if (!File::Exists(outFile))
        {
            MessageBox::Show("Файл результатов не найден: " + outFile);
            return;
        }

        cli::array<String^>^ lines = File::ReadAllLines(outFile);
        List<double>^ Xp = gcnew List<double>();
        List<double>^ P = gcnew List<double>();

        bool readXp = false;
        bool readP = false;

        for each(String ^ line in lines)
        {
            String^ trimmed = line->Trim();

            if (trimmed == "Xp")
            {
                readXp = true;
                readP = false;
                continue;
            }
            else if (trimmed == "P")
            {
                readP = true;
                readXp = false;
                continue;
            }
            else if (trimmed == "" || trimmed->Contains("="))
            {
                readXp = false;
                readP = false;
                continue;
            }

            if (readXp || readP)
            {
                cli::array<String^>^ values = trimmed->Split(
                    gcnew cli::array<wchar_t>{' ', '\t', ',', ';'},
                    StringSplitOptions::RemoveEmptyEntries);

                for each(String ^ value in values)
                {
                    double num;
                    if (Double::TryParse(value, num))
                    {
                        if (readXp) Xp->Add(num);
                        else if (readP) P->Add(num);
                    }
                }
            }
        }

        if (Xp->Count > 0 && P->Count > 0 && Xp->Count == P->Count)
        {
            for (int i = 0; i < P->Count - 1; i++)
            {
                for (int j = i + 1; j < P->Count; j++)
                {
                    if (P[i] > P[j])
                    {
                        double tempP = P[i];
                        double tempX = Xp[i];
                        P[i] = P[j];
                        Xp[i] = Xp[j];
                        P[j] = tempP;
                        Xp[j] = tempX;
                    }
                }
            }

            Series^ distSeries = gcnew Series("Функция распределения");
            distSeries->ChartType = SeriesChartType::Line;
            distSeries->Color = Color::Blue;
            distSeries->BorderWidth = 3;
            distSeries->MarkerSize = 6;
            distSeries->MarkerStyle = MarkerStyle::Circle;
            distSeries->MarkerColor = Color::Red;

            for (int i = 0; i < P->Count; i++)
            {
                distSeries->Points->AddXY(P[i], Xp[i]);
            }

            chart->Series->Add(distSeries);

            chart->ChartAreas[0]->AxisX->Title = "Вероятность P";
            chart->ChartAreas[0]->AxisY->Title = "Квантиль Xp";
            chart->ChartAreas[0]->AxisX->Minimum = 0;
            chart->ChartAreas[0]->AxisX->Maximum = 1;

            UpdateStatus("Построен график функции распределения: " + method);
        }
        else
        {
            MessageBox::Show("Не удалось найти данные Xp и P в файле результатов");
        }
    }


private:
    void UpdateStatus(String^ message)
    {
        statusStrip->Items->Clear();
        statusStrip->Items->Add(message);
    }

    void ShowDataChart(String^ dataFile, String^ title)
    {
        try
        {
            StreamReader^ reader = gcnew StreamReader(dataFile);
            Series^ series = gcnew Series();
            series->ChartType = SeriesChartType::Line;
            series->Color = Color::Blue;
            series->BorderWidth = 2;
            series->Name = "Data";

            String^ line;
            while ((line = reader->ReadLine()) != nullptr)
            {
                if (String::IsNullOrEmpty(line) || line->StartsWith("#"))
                    continue;

                Collections::Generic::List<String^>^ parts = gcnew Collections::Generic::List<String^>();
                int start = 0;
                for (int i = 0; i <= line->Length; i++)
                {
                    if (i == line->Length || line[i] == ' ' || line[i] == '\t')
                    {
                        if (start < i)
                        {
                            String^ part = line->Substring(start, i - start);
                            parts->Add(part);
                        }
                        start = i + 1;
                    }
                }

                if (parts->Count >= 2)
                {
                    double x, y;
                    if (Double::TryParse(parts[0], x) && Double::TryParse(parts[1], y))
                    {
                        series->Points->AddXY(x, y);
                    }
                }
            }
            reader->Close();

            chart->Series->Add(series);
            chart->Titles->Clear();
            chart->Titles->Add(title + " - " + currentMethod);
            UpdateStatus("График построен: " + currentMethod);
        }
        catch (Exception^ ex)
        {
            MessageBox::Show("Ошибка построения графика: " + ex->Message, "Ошибка");
        }
    }

    void ShowQQChart(String^ dataFile)
    {
        try
        {
            StreamReader^ reader = gcnew StreamReader(dataFile);

            Series^ pointsSeries = gcnew Series();
            pointsSeries->ChartType = SeriesChartType::Point;
            pointsSeries->Color = Color::Blue;
            pointsSeries->MarkerSize = 6;
            pointsSeries->Name = "QQ Points";

            Series^ lineSeries = gcnew Series();
            lineSeries->ChartType = SeriesChartType::Line;
            lineSeries->Color = Color::Red;
            lineSeries->BorderWidth = 2;
            lineSeries->Name = "Reference";

            double minX = Double::MaxValue;
            double maxX = Double::MinValue;
            double minY = Double::MaxValue;
            double maxY = Double::MinValue;

            String^ line;
            while ((line = reader->ReadLine()) != nullptr)
            {
                if (String::IsNullOrEmpty(line) || line->StartsWith("#"))
                    continue;

                Collections::Generic::List<String^>^ parts = gcnew Collections::Generic::List<String^>();
                int start = 0;
                for (int i = 0; i <= line->Length; i++)
                {
                    if (i == line->Length || line[i] == ' ' || line[i] == '\t')
                    {
                        if (start < i)
                        {
                            String^ part = line->Substring(start, i - start);
                            parts->Add(part);
                        }
                        start = i + 1;
                    }
                }

                if (parts->Count >= 2)
                {
                    double x, y;
                    if (Double::TryParse(parts[0], x) && Double::TryParse(parts[1], y))
                    {
                        pointsSeries->Points->AddXY(x, y);
                        minX = Math::Min(minX, x);
                        maxX = Math::Max(maxX, x);
                        minY = Math::Min(minY, y);
                        maxY = Math::Max(maxY, y);
                    }
                }
            }
            reader->Close();

            if (pointsSeries->Points->Count > 0)
            {
                double minVal = Math::Min(minX, minY);
                double maxVal = Math::Max(maxX, maxY);
                lineSeries->Points->AddXY(minVal, minVal);
                lineSeries->Points->AddXY(maxVal, maxVal);
            }

            chart->Series->Add(pointsSeries);
            chart->Series->Add(lineSeries);
            chart->Titles->Clear();
            chart->Titles->Add("Q-Q график - " + currentMethod);
            UpdateStatus("Q-Q график построен: " + currentMethod);
        }
        catch (Exception^ ex)
        {
            MessageBox::Show("Ошибка построения Q-Q графика: " + ex->Message, "Ошибка");
        }
    }
};
