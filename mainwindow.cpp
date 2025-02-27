#include "mainwindow.h"
#include "components.h"
#include "ui_mainwindow.h"
#include <fftw3.h>
#include <QVector>

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
    laser.frequencyResolution = ui->lineEdit_res->text().toDouble();

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
        double freqMin = spectrumData.frequency.first();
        double freqMax = spectrumData.frequency.last();
        double dataRange = freqMax - freqMin;
        // Если диапазон данных слишком узкий, задаём минимальный видимый диапазон оси X:
        double desiredRange = (dataRange > 0) ? std::max(dataRange, dataRange / 0.4) : 1.0;
        double center = (freqMin + freqMax) / 2.0;

        ui->pulse_plot->xAxis->setRange(center - desiredRange / 2, center + desiredRange / 2);
        ui->pulse_plot->yAxis->setRange(0, *std::max_element(spectrumData.intensity.begin(), spectrumData.intensity.end()));
    }


    ui->pulse_plot->replot();
}

void MainWindow::plotTimeDomain(const Laser &laser) {
    Components components;
    SpectrumData spectrumData = components.get_spectrum(laser);
    TimeDomainData timeData = components.get_time_domain(spectrumData, laser);

    ui->pulse_plot->clearGraphs();
    ui->pulse_plot->addGraph();
    ui->pulse_plot->graph(0)->setData(timeData.time, timeData.intensity);

    ui->pulse_plot->xAxis->setLabel("Время (с)");
    ui->pulse_plot->yAxis->setLabel("Интенсивность");

    if (!timeData.time.isEmpty() && !timeData.intensity.isEmpty()) {
        ui->pulse_plot->xAxis->setRange(timeData.time.first(), timeData.time.last());
        ui->pulse_plot->yAxis->setRange(0, *std::max_element(timeData.intensity.begin(), timeData.intensity.end()));
    }

    ui->pulse_plot->replot();
}


