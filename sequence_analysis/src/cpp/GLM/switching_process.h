#ifndef SWITCHING_PROCESS_H
#define SWITCHING_PROCESS_H

#include <sstream>
#include "discrete_parametric.h"


const double SWITCHING_NOISE_PROBABILITY = 0.05; // perturbation des fonctions de masse  des observations

class Switching_process {

  friend std::ostream& operator<<(std::ostream &os, const Switching_process &sw_process)
  { return sw_process.ascii_print(os, 0, false, false);}

 public:
  int nb_state; // nombre d'états
  int nb_value; // nombre de valeurs
  Discrete_parametric **observation; //lois d'observation

  Switching_process(int inb_state = 0, int inb_value = 0);
  Switching_process(int inb_state, Discrete_parametric **pobservation);
  Switching_process(const Switching_process &process, char manip = 'c', int state = -1);
  ~Switching_process();

  Switching_process& operator=(const Switching_process &process);

  void copy(const Switching_process  &process);
  void add_state(const Switching_process &process, int state);
  void remove();

  std::ostream& ascii_print(std::ostream &os, Histogram **empirical_observation = 0 ,
			    bool exhaustive = true, bool file_flag = true) const;
  std::ostream& spreadsheet_print(std::ostream &os , Histogram **empirical_observation = 0) const;

  void nb_value_computation();
  int nb_parameter_computation() const;
  
  void init();

  // private:
  // protected:
};

#endif
