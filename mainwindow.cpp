#include "mainwindow.h"
#include "components.h"
#include "ui_mainwindow.h"
#include <fftw3.h>
#include <QVector>
#include <QDebug>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->pulse_plot->setInteraction(QCP::iRangeZoom, true);
    ui->pulse_plot->setInteraction(QCP::iRangeDrag, true);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_pushButton_Start_clicked()
{
    laser.centralWavelength = ui->lineEdit_centralWavelength->text().toDouble();
    laser.phase = ui->lineEdit_phase->text().toDouble();
    laser.pulseDuration = ui->lineEdit_pulseDuration->text().toDouble();
    laser.averageCountPhotons = ui->lineEdit_averageCountPh->text().toDouble();
    laser.numberPoints = ui->lineEdit_N->text().toDouble();

    // Проверяем, какой QRadioButton выбран
    if (ui->radioButton_spec->isChecked()) {
        // Построение спектра
        plotGraph(laser);
    }
    else if (ui->radioButton_time->isChecked()) {
        // Построение временной области на основе спектра
        plotTimeDomain(laser);
    }
}

void MainWindow::plotGraph(const Laser &laser) {
    Components components;
    SpectrumData spectrumData = components.get_spectrum(laser);

    // Очищаем предыдущие графики
    ui->pulse_plot->clearGraphs();
    ui->pulse_plot->addGraph();

    // Устанавливаем данные для графика
    ui->pulse_plot->graph(0)->setData(spectrumData.frequency, spectrumData.intensity);

    ui->pulse_plot->xAxis->setLabel("Частота (Гц)");
    ui->pulse_plot->yAxis->setLabel("Амплитуда");

    // Установка диапазонов осей
    if (!spectrumData.frequency.isEmpty() && !spectrumData.intensity.isEmpty()) {
        double freq_min = spectrumData.frequency.first();
        double freq_max = spectrumData.frequency.last();
        double data_range = freq_max - freq_min;
        // Если диапазон данных слишком узкий, задаём минимальный видимый диапазон оси X:
        double desired_range = (data_range > 0) ? std::max(data_range, data_range / 0.4) : 1.0;
        double center = (freq_min + freq_max) / 2.0;

        ui->pulse_plot->xAxis->setRange(center - desired_range / 2, center + desired_range / 2);
        ui->pulse_plot->yAxis->setRange(0, *std::max_element(spectrumData.intensity.begin(), spectrumData.intensity.end()));
    }

    ui->pulse_plot->replot();
}

TimeDomainData generateCompositePulse(const TimeDomainData &singlePulse, int numPulses, double dt)
{
    // Количество точек в одном импульсе
    int N_single = singlePulse.time.size();
    if (N_single == 0 || numPulses < 1)
        return singlePulse; // защита от некорректных данных

    // Рассчитаем длительность одного импульса
    double T_pulse = singlePulse.time.last() - singlePulse.time.first();
    // Зададим временной интервал между началом соседних импульсов.
    // Здесь выбран коэффициент 1.5 (можно изменить по необходимости).
    double T_sep = T_pulse * 1.5;
    int shiftSamples = static_cast<int>(std::round(T_sep / dt));

    // Общая длина результирующего сигнала
    int N_composite = N_single + (numPulses - 1) * shiftSamples;

    TimeDomainData composite;
    composite.time.resize(N_composite);
    composite.intensity.resize(N_composite, 0.0);

    // Заполним временную ось результирующего сигнала
    double t0 = singlePulse.time.first();
    for (int i = 0; i < N_composite; ++i) {
        composite.time[i] = t0 + i * dt;
    }

    // Суммируем копии одиночного импульса, сдвигая их на shiftSamples
    for (int p = 0; p < numPulses; ++p) {
        int offset = p * shiftSamples;
        for (int j = 0; j < N_single; ++j) {
            if (offset + j < N_composite)
                composite.intensity[offset + j] += singlePulse.intensity[j];
        }
    }
    return composite;
}

void MainWindow::plotTimeDomain(const Laser &laser)
{
    Components components;
    // 1) Считаем спектр и преобразуем его в временную область
    SpectrumData spectrumData = components.get_spectrum(laser);
    TimeDomainData singlePulse = components.get_time_domain(spectrumData, laser);

    // 2) Читаем число импульсов из lineEdit_num_pulse (по умолчанию должно быть 1 или больше)
    int numPulses = ui->lineEdit_num_pulse->text().toInt();
    // Если число импульсов больше 1, формируем композицию
    TimeDomainData timeData;
    if(numPulses > 1 && singlePulse.time.size() >= 2) {
        double dt = singlePulse.time[1] - singlePulse.time[0];
        timeData = generateCompositePulse(singlePulse, numPulses, dt);
    }
    else {
        timeData = singlePulse;
    }

    // 3) Строим график во временной области
    ui->pulse_plot->clearGraphs();
    ui->pulse_plot->addGraph();
    ui->pulse_plot->graph(0)->setData(timeData.time, timeData.intensity);
    ui->pulse_plot->xAxis->setLabel("Время (с)");
    ui->pulse_plot->yAxis->setLabel("Интенсивность (Вт)");

    if(!timeData.time.isEmpty() && !timeData.intensity.isEmpty()) {
        double tMin = timeData.time.first();
        double tMax = timeData.time.last();
        double iMax = *std::max_element(timeData.intensity.begin(), timeData.intensity.end());
        ui->pulse_plot->xAxis->setRange(tMin, tMax);
        ui->pulse_plot->yAxis->setRange(0, iMax);
    }
    ui->pulse_plot->replot();

    // 4) Подсчёт энергии одного импульса.
    // Здесь энергия считается как интеграл I(t) dt для всей композиции,
    // делённый на число импульсов.
    double totalEnergy = 0.0;
    int N_time = timeData.time.size();
    for (int i = 0; i < N_time - 1; ++i) {
        double dt = timeData.time[i+1] - timeData.time[i];
        double I_avg = 0.5 * (timeData.intensity[i] + timeData.intensity[i+1]);
        totalEnergy += I_avg * dt;
    }
    double pulseEnergy = totalEnergy / numPulses;
    qDebug() << "Energy of one pulse:" << pulseEnergy << "J";
}
