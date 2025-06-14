#include "mainwindow.h"
#include "components.h"
#include "ui_mainwindow.h"
#include <QVector>
#include <QDebug>
#include <QFile>
#include <QTextStream>
#include <QDir>
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
    laser.central_wavelength = ui->lineEdit_centralWavelength->text().toDouble();
    laser.pulse_duration = ui->lineEdit_pulseDuration->text().toDouble();
    laser.avg_count_photons = ui->lineEdit_averageCountPh->text().toDouble();
    laser.number_points = ui->lineEdit_N->text().toDouble();
    laser.rep_rate = ui->lineEdit_repRate->text().toDouble(); // частота в МГц

    quantumChannel.channel_attenuation = ui->lineEdit_att->text().toDouble();
    quantumChannel.channel_length = ui->lineEdit_length->text().toDouble();
    quantumChannel.chromatic_dispersion = ui->lineEdit_crom_disp->text().toDouble();

    detector.quantum_efficiency = ui->lineEdit_QE->text().toDouble();
    detector.dead_time = ui->lineEdit_dead_time->text().toDouble();
    detector.time_slot = ui->lineEdit_time_slot->text().toDouble();

    int num_repeat_experiment = ui->lineEdit_repeat->text().toInt();
    double step_PD = ui->lineEdit_stepPD->text().toDouble();
    double step_RR = ui->lineEdit_stepRR->text().toDouble();

    // Выбор режима построения графика
    if (ui->radioButton_spec->isChecked()) {
        // Построение спектра
        plotGraph(laser);
    }
    else if (ui->radioButton_time->isChecked()) {
        // Построение временной области с генерацией серии импульсов
        plotTimeDomain(laser);
    }
    else if (ui->radioButton_gen_key->isChecked() && num_repeat_experiment == 1){
        plotGenKeys(laser, detector);
    }
    else if (ui->radioButton_gen_key->isChecked() && num_repeat_experiment > 1){
        for (int i = 0; i < num_repeat_experiment; i++) {
            laser.pulse_duration = ui->lineEdit_pulseDuration->text().toDouble() - i * step_PD;
            laser.rep_rate = ui->lineEdit_repRate->text().toDouble() + i * step_RR; // частота в МГц
            plotGenKeys(laser, detector);
        }
        qDebug() << "End";
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
    // Получаем спектр и преобразуем его во временную область
    SpectrumData spectrumData = components.get_spectrum(laser);
    TimeDomainData singlePulse = components.spectrum_to_time_domain(spectrumData, laser, quantumChannel);

    int numPulses = ui->lineEdit_num_pulse->text().toInt();
    if (numPulses < 1) {
        numPulses = 1;
    }

    TimeDomainData timeData;
    if (numPulses > 1 && singlePulse.time.size() >= 2) {
        double dt = singlePulse.time[1] - singlePulse.time[0];
        // Вызов метода для построения графика во временной обл
        timeData = components.gen_composite_pulse(singlePulse, laser, numPulses, dt, quantumChannel);
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

    // Строим график во временной области
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

    // Подсчёт энергии одного импульса (интеграл I(t) dt, делённый на число импульсов)
    double totalEnergy = 0.0;
    int N_time = timeData.time.size();
    for (int i = 0; i < N_time - 1; ++i) {
        double dt = timeData.time[i+1] - timeData.time[i];
        double I_avg = 0.5 * (timeData.intensity[i] + timeData.intensity[i+1]);
        totalEnergy += I_avg * dt;
    }
//    double pulseEnergy = totalEnergy / numPulses;
//    qDebug() << "Energy of one pulse:" << pulseEnergy << "J";
}

void MainWindow::plotGenKeys(const Laser &laser, const Photodetector &detector)
{
    Components components;
    // Получаем спектр и преобразуем его во временную область
    SpectrumData spectrumData = components.get_spectrum(laser);
    TimeDomainData singlePulse = components.spectrum_to_time_domain(spectrumData, laser, quantumChannel);
    std::vector<unsigned int> photon_counts; // стартовый вектор импульсов
    std::vector<std::vector<double>> matrix_pulses; // матрица импульсов с фотонами

    unsigned int num_pulses = ui->lineEdit_num_pulse->text().toInt();
    if (num_pulses > 0)
    {
        for (unsigned int i = 0; i < num_pulses; i++) {
            unsigned int ph = components.get_photons(laser, quantumChannel);
            photon_counts.push_back(ph); // заполнение импульсами
        }
    }
    //qDebug() << "Vector:" << photon_counts;

    // время импульса
    std::vector<double> time_pulse;
    std::vector<double> intensity_pulse;
    for (int i = 0; i < singlePulse.time.size(); i++) {
        if(singlePulse.intensity[i] > 2e-12){
            time_pulse.push_back(singlePulse.time[i]);
            intensity_pulse.push_back(singlePulse.intensity[i]);
        }
    }
//    qDebug() << "Time (vector):" << time_pulse;
    double time = 0;
    for (unsigned int i = 0; i < time_pulse.size() ;i++ ) {
        time = abs(time_pulse[i]) + time;
    }
//    qDebug() << "Time:" << time;

    std::vector<double> time_slots;
    components.get_time_slot(num_pulses, laser, time_slots, detector);
//    qDebug() << "Time slot:" << time_slots;

    // временные метки
    components.gen_ph_timelabel(num_pulses, photon_counts, matrix_pulses, detector, time_slots, time);
    std::vector<double> time_reg;
    int num_reg_pulses = components.reg_pulses(matrix_pulses, detector, time_slots, time_reg);

    if(time_reg.size() <= 1)
    {
//        qDebug() << "Не было зарегистрированных импульсов";
    }
    else{
//        qDebug() << "Кол-во зарег. импульсов: " << num_reg_pulses;
        ui->lineEdit_count_pulses->clear();
        ui->lineEdit_percent_pulses->clear();
        ui->lineEdit_count_pulses->setText(QString::number(num_reg_pulses)); // число принятых импульсов
        double result = (num_reg_pulses * 1.0) / num_pulses;
        ui->lineEdit_percent_pulses->setText(QString::number(result, 'f', 2)); // отношение принятых имп к отправ.

//        qDebug() << "Vector:" << time_reg;
        QVector<double> value_pulse;
        value_pulse.resize(num_reg_pulses);
        QVector<double> time_pulse;
        time_pulse.resize(num_reg_pulses);

        for (int i = 0; i < num_reg_pulses; i++ ) {
            value_pulse.fill(1);
            time_pulse[i] = time_reg[i + 1];
        }
//        qDebug() << value_pulse;
//        qDebug() << time_pulse;

        // вывод в файл
        save_results_to_file(num_reg_pulses, num_pulses,
                             laser.pulse_duration,
                             laser.rep_rate,
                             detector.quantum_efficiency,
                             detector.dead_time,
                             detector.time_slot);

        ui->pulse_plot->clearGraphs();

        QCPGraph *stemPlot = ui->pulse_plot->addGraph();
        stemPlot->setLineStyle(QCPGraph::lsImpulse);

        stemPlot->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 5));
        stemPlot->setData(time_pulse, value_pulse);

        ui->pulse_plot->xAxis->setLabel("Время (с)");
        ui->pulse_plot->yAxis->setLabel("Число фотонов");
        if (!time_pulse.isEmpty() && !value_pulse.isEmpty()) {
            double tMin = time_pulse.first();
            double tMax = time_pulse.last();
            double iMax = *std::max_element(value_pulse.begin(), value_pulse.end());
            ui->pulse_plot->xAxis->setRange(tMin + tMin * 0.2, tMax + tMax * 0.2);
            ui->pulse_plot->yAxis->setRange(0, iMax + 0.3);
        }
        ui->pulse_plot->replot();
    }
}

void MainWindow::save_results_to_file(int num_reg_pulses, int total_sent_pulses, double pulse_duration,
                                      double repRate, double QE, double dead_time, double time_slot) {
    double result = (double)num_reg_pulses / total_sent_pulses;

    QString filename = QDir::currentPath() + "/simulation_results.csv";
    QFile file(filename);
    bool fileExists = file.exists();

    if (file.open(QIODevice::Append | QIODevice::Text)) {
        QTextStream out(&file);

        // Заголовки только если файл создаётся впервые
        if (!fileExists) {
            out << "pulse_duration_ns,repRate_MHz,quantum_eff,dead_time_ns,time_slot_ns,num_registered,total_sent,result\n";
        }

        out << pulse_duration << "," << repRate << "," << QE << "," << dead_time << "," << time_slot << ","
            << num_reg_pulses << "," << total_sent_pulses << "," << result << "\n";

        file.close();
//        qDebug() << "Результаты сохранены в" << filename;
    } else {
        qDebug() << "Ошибка открытия файла для записи результатов!";
    }
}
