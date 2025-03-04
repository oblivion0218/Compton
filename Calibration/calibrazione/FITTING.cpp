#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"

using namespace std;

// Funzione per leggere i dati di un file con posizione picchi e corrispettiva altezza
void read_fit_results(const string& filename, vector<int>& peaks, vector<double>& centers, vector<double>& heights) {
    ifstream file(filename);
    string line;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Salta righe vuote o commenti
        }
        istringstream iss(line);
        int peak;
        double center, height;
        char comma;
        iss >> peak >> comma >> center >> comma >> height;
        peaks.push_back(peak);
        centers.push_back(center);
        heights.push_back(height);
    }
    file.close();
}

// funzione rubata ad ANDRE e riadattata
vector<int> read_histogram_data(const string& filename) {
    ifstream file(filename);
    
    vector<int> histogram_data;
    string line;
    bool found_data_section = false;

    // Cerca la sezione $DATA
    while (getline(file, line)) {
        if (line.find("$DATA:") != string::npos) {
            found_data_section = true;
            break;
        }
    }

     // Salta la riga con l'intervallo "0 2047"
    if (!std::getline(file, line)) {
        throw std::runtime_error("Errore nella lettura della riga dell'intervallo");
    }

    // Legge i valori dei bin dell'istogramma
    while (std::getline(file, line)) {
        try {
            int bin_value = std::stoi(line);
            histogram_data.push_back(bin_value);
        } catch (const std::invalid_argument&) {
            break;
        }
    }

    return histogram_data;
}



int main() {

    string name = "44-Ti";
    cout << "Lettura dati da " << name << "..." << endl;

    // Lettura dei dati dall'istogramma
    vector<int> histogram_data = read_histogram_data(name + ".Spe");

    // Creazione dell'oggetto TH1F
    int n_bins = histogram_data.size();
    TH1F* hist = new TH1F("hist", "", n_bins, 0, n_bins);
    for (size_t i = 0; i < histogram_data.size(); ++i) {
        hist->SetBinContent(i + 1, histogram_data[i]);
    }

    gStyle->SetOptStat(0); // Disabilita il box con i conteggi e deviazione standard
    hist->SetStats(0);

    // Lettura dei risultati del fit
    vector<int> peaks;
    vector<double> centers;
    vector<double> heights;
    read_fit_results(name + "_fit_results.txt", peaks, centers, heights);

    vector<vector<double>> fit_results;

    for (size_t i = 0; i < peaks.size(); ++i) {

        cout << "Eseguo il fit per il picco alla posizione " << peaks[i] << "..." << endl;

        // Definire la finestra di fit intorno al picco
        double extreme = 20; // Ampiezza di fit intorno al picco
        double fit_min = centers[i] - extreme;
        double fit_max = centers[i] + extreme;

        cout <<" \n min : " << fit_min << "\n max: " << fit_max<<endl;

        // Definizione della funzione di fit
        TF1* f_gauss = new TF1("gauss+fondo", "gaus(0) + pol2(3)", fit_min, fit_max);
        f_gauss->SetParameters(heights[i] ,centers[i], 10, heights[i+1], centers[i+1], 10); // Parametri iniziali del fit
        f_gauss->SetParLimits(0, heights[i]-100 , heights[i]+100);
        f_gauss->SetParLimits(1, centers[i]-100 , centers[i]+100);


        // Eseguire il fit sull'istogramma
        hist->Fit(f_gauss, "RQ");

        // Estrarre i parametri del fit e gli errori
        double true_height = f_gauss->GetParameter(0);
        double height_error = f_gauss->GetParError(0);
        double true_center = f_gauss->GetParameter(1);
        double center_error = f_gauss->GetParError(1);
        double width = f_gauss->GetParameter(2);
        double width_error = f_gauss->GetParError(2);

        double chi2 = f_gauss->GetChisquare(); // Chi-quadro del fit
        double ndf = f_gauss->GetNDF();        // Numero di gradi di libertà
        double chi2_reduced = (ndf > 0) ? chi2 / ndf : 0; // Calcola il chi-quadro ridotto

       
        // Stampare il grafico del fit con istogramma limitato all'intervallo di fit
        TCanvas* canvas = new TCanvas("c1", "c1", 800, 600);
        canvas->SetTitle((name + "_" + to_string(centers[i])).c_str());
        hist->GetXaxis()->SetRangeUser(fit_min - 100, fit_max + 100); // Limitare l'intervallo dell'asse x
        hist->Draw(); // Disegnare l'istogramma
        f_gauss->Draw("SAME"); // Sovrapporre il fit

        // Aggiungere i parametri del fit al grafico
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.03);
        latex.DrawLatex(0.15, 0.85, Form("Peak: %d", peaks[i]));
        latex.DrawLatex(0.15, 0.80, Form("#chi^{2}/ndf: %.2f", chi2_reduced));
        latex.DrawLatex(0.15, 0.75, Form("Center1: %.2f +/- %.2f", true_center, center_error));
        latex.DrawLatex(0.15, 0.70, Form("Height1: %.2f +/- %.2f", true_height, height_error));
        latex.DrawLatex(0.15, 0.65, Form("Width1: %.2f +/- %.2f", width, width_error));

        // Salvare l'immagine del fit
        canvas->SaveAs((name + "_"  + to_string(i+1) + "_fit.png").c_str());

        delete canvas;
        delete f_gauss;
    }

    delete hist;
    return 0;
}
