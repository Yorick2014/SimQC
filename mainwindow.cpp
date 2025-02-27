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

void MainWindow::plotTimeDomain(const Laser &laser)
{
    Components components;
    // 1) Считаем спектр
    SpectrumData spectrumData = components.get_spectrum(laser);

    // 2) Преобразуем в временную область (с нормировкой на photon count)
    TimeDomainData timeData = components.get_time_domain(spectrumData, laser);

    // 3) Строим график
    ui->pulse_plot->clearGraphs();
    ui->pulse_plot->addGraph();
    ui->pulse_plot->graph(0)->setData(timeData.time, timeData.intensity);
    ui->pulse_plot->xAxis->setLabel("Время (с)");
    ui->pulse_plot->yAxis->setLabel("Интенсивность (Вт)");

    if (!timeData.time.isEmpty() && !timeData.intensity.isEmpty()) {
        double time_min = timeData.time.first();
        double time_max = timeData.time.last();
        double intensity_max = *std::max_element(timeData.intensity.begin(), timeData.intensity.end());
        ui->pulse_plot->xAxis->setRange(time_min, time_max);
        ui->pulse_plot->yAxis->setRange(0, intensity_max);
    }

    ui->pulse_plot->replot();

    // Подсчёт энергии импульса: интеграл от I(t) по времени
    double pulse_energy = 0.0;
    int N_time = timeData.time.size();
    for (int i = 0; i < N_time - 1; ++i) {
        double dt = timeData.time[i+1] - timeData.time[i];
        // Метод трапеций для аппроксимации интеграла
        double intensity_avg = 0.5 * (timeData.intensity[i] + timeData.intensity[i+1]);
        pulse_energy += intensity_avg * dt;
    }

    // Вывод энергии импульса (в Джоулях) в терминал
    qDebug() << "Pulse energy:" << pulse_energy << "J";
}
