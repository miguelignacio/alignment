#ifndef _plot_h
#define _plot_h

// ROOT
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TStyle.h"


// event
#include "AlignEvent.h"


// Globals
//=========

/// language in plots
extern bool gGerman;

// Method prototypes
//===================

char * strdup_new(const char * text);

// Standard option settings
void setopt(TStyle * style);
void setopt(TCanvas * canvas);
void setopt(TVirtualPad * pad);
void setopt(TH1 * histo);
void setopt(TLegend *leg);
void setopt(TGraph *gr);
void setoptlarge(TH1 * histo);

const char * GetUnit(const char * hname);
const char * GetXTitle(const TH1F * histo);
const char * GetYTitle(const TH1F * histo);
void settitle(TH1F * histo);

#endif
