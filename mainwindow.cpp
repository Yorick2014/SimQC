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

    plotGraph(laser);
}

void MainWindow::plotGraph(const Laser &laser)
{
    Components components;
    std::vector<double> frequency, spectrum;

    // Вызываем функцию расчёта спектра
    components.spectrum(laser, frequency, spectrum);

    // Преобразуем std::vector в QVector для графика
    QVector<double> qFrequency = QVector<double>(frequency.begin(), frequency.end());
    QVector<double> qSpectrum = QVector<double>(spectrum.begin(), spectrum.end());

    // Удаляем предыдущие графики
    ui->pulse_plot->clearGraphs();
    ui->pulse_plot->addGraph();

    // Добавляем новый график
    ui->pulse_plot->graph(0)->setData(qFrequency, qSpectrum);

    ui->pulse_plot->xAxis->setLabel("Частота (Гц)");
    ui->pulse_plot->yAxis->setLabel("Амплитуда");

    // Установка диапазонов осей
    if (!qFrequency.isEmpty() && !qSpectrum.isEmpty()) {
        ui->pulse_plot->xAxis->setRange(0, qFrequency.last());
        ui->pulse_plot->yAxis->setRange(0, *std::max_element(qSpectrum.begin(), qSpectrum.end()));
    }

    ui->pulse_plot->replot();
}

