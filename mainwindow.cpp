#include "mainwindow.h"
#include "components.h"
#include "ui_mainwindow.h"
#include <QVector>
#include <QDebug>
#include <cmath>
#include <random>
#include <iomanip>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow)
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
    quantumChannel.isAtt = false;
    quantumChannel.isCromDisp = false;
    if (ui->checkBox_att->isChecked()){
        quantumChannel.isAtt = true;
    }
    if (ui->checkBox_disp->isChecked()){
        quantumChannel.isCromDisp = true;
    }

    // Считываем параметры лазера
    laser.centralWavelength = ui->lineEdit_centralWavelength->text().toDouble();
    laser.phase = ui->lineEdit_phase->text().toDouble();
    laser.pulseDuration = ui->lineEdit_pulseDuration->text().toDouble();
    laser.averageCountPhotons = ui->lineEdit_averageCountPh->text().toDouble();
    laser.numberPoints = ui->lineEdit_N->text().toDouble();
    laser.repRate = ui->lineEdit_repRate->text().toDouble(); // частота в МГц

    quantumChannel.channelAttenuation = ui->lineEdit_att->text().toDouble();
    quantumChannel.channelLength = ui->lineEdit_length->text().toDouble();
    quantumChannel.chromaticDispersion = ui->lineEdit_crom_disp->text().toDouble();

    // Выбор режима построения графика
    if (ui->radioButton_spec->isChecked()) {
        // Построение спектра
        plotGraph(laser);
    }
    else if (ui->radioButton_time->isChecked()) {
        // Построение временной области с генерацией серии импульсов
        plotTimeDomain(laser);
    }
}

void MainWindow::plotGraph(const Laser &laser) {
    Components components;
    SpectrumData spectrumData = components.get_spectrum(laser);

    ui->pulse_plot->clearGraphs();
    ui->pulse_plot->addGraph();
    ui->pulse_plot->graph(0)->setData(spectrumData.frequency, spectrumData.intensity);

    ui->pulse_plot->xAxis->setLabel("Частота (Гц)");
    ui->pulse_plot->yAxis->setLabel("Амплитуда");

    if (!spectrumData.frequency.isEmpty() && !spectrumData.intensity.isEmpty()) {
        double freq_min = spectrumData.frequency.first();
        double freq_max = spectrumData.frequency.last();
        double data_range = freq_max - freq_min;
        double desired_range = (data_range > 0) ? std::max(data_range, data_range / 0.4) : 1.0;
        double center = (freq_min + freq_max) / 2.0;

        ui->pulse_plot->xAxis->setRange(center - desired_range / 2, center + desired_range / 2);
        ui->pulse_plot->yAxis->setRange(0, *std::max_element(spectrumData.intensity.begin(), spectrumData.intensity.end()));
    }
    ui->pulse_plot->replot();
}

void MainWindow::plotTimeDomain(const Laser &laser)
{
    Components components;
    // 1) Получаем спектр и преобразуем его во временную область (одиночный импульс)
    SpectrumData spectrumData = components.get_spectrum(laser);
    TimeDomainData singlePulse = components.get_time_domain(spectrumData, laser, quantumChannel);

    // 2) Читаем число импульсов (из lineEdit_num_pulse)
    int numPulses = ui->lineEdit_num_pulse->text().toInt();
    if (numPulses < 1) {
        numPulses = 1;
    }

    TimeDomainData timeData;
    if (numPulses > 1 && singlePulse.time.size() >= 2) {
        double dt = singlePulse.time[1] - singlePulse.time[0];
        // Вызов метода для построения графика во временной обл
        timeData = components.generateCompositePulse(singlePulse, laser, numPulses, dt, quantumChannel);
    }
    else {
        timeData = singlePulse;
    }

    // Найдем первый индекс, где интенсивность становится больше 0
    int idx = -1;
    for (int i = 0; i < timeData.intensity.size(); ++i) {
        if (timeData.intensity[i] > 2e-14) {
            idx = i;
            break;
        }
    }

    if (idx != -1) {
        // Получаем время, соответствующее первому ненулевому значению интенсивности
        double shiftTime = timeData.time[idx];

        // Смещаем все значения времени так, чтобы время в этом индексе стало 0
        for (int i = 0; i < timeData.time.size(); ++i) {
            timeData.time[i] -= shiftTime;
        }
    }

    // 3) Строим график во временной области
    ui->pulse_plot->clearGraphs();
    ui->pulse_plot->addGraph();
    ui->pulse_plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 3));
    ui->pulse_plot->graph(0)->setData(timeData.time, timeData.intensity);
    ui->pulse_plot->xAxis->setLabel("Время (с)");
    ui->pulse_plot->yAxis->setLabel("Интенсивность (Вт)");

    if (!timeData.time.isEmpty() && !timeData.intensity.isEmpty()) {
        double tMin = timeData.time.first();
        double tMax = timeData.time.last();
        double iMax = *std::max_element(timeData.intensity.begin(), timeData.intensity.end());
        ui->pulse_plot->xAxis->setRange(tMin, tMax);
        ui->pulse_plot->yAxis->setRange(0, iMax);
    }
    ui->pulse_plot->replot();

    // 4) Подсчёт энергии одного импульса (интеграл I(t) dt, делённый на число импульсов)
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
