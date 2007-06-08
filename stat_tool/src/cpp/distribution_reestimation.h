/*--------------------------------------------------------------*
 *
 *  Reestimation des probabilites d'une loi discrete.
 *
 *  arguments : nombre de valeurs, pointeurs sur les quantites
 *              de reestimation et sur les probabilites de la loi discrete,
 *              probabilite minimum, flag pour reestimer des probabilites nulles.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void reestimation(int nb_value , Type *reestim , double *pmass ,
                  double min_probability , bool null_probability)

{
  register int i;
  int nb_correction;
  double sum , norm;


  sum = 0.;
  for (i = 0;i < nb_value;i++) {
    sum += *reestim++;
  }
  reestim -= nb_value;

  if (sum > 0.) {
    norm = 0.;
    nb_correction = 0;

    for (i = 0;i < nb_value;i++) {
      if (*pmass > 0.) {
        if (*reestim > min_probability * sum) {
          *pmass = *reestim / sum;
          norm += *pmass;
        }

        else {
          if ((*reestim > 0.) || (!null_probability)) {
            nb_correction++;
            *pmass = min_probability;
          }
          else {
            *pmass = 0.;
          }
        }
      }

      reestim++;
      pmass++;
    }

    if (nb_correction > 0) {
      reestim -= nb_value;
      pmass -= nb_value;

      for (i = 0;i < nb_value;i++) {
        if ((*pmass > 0.) && (*reestim > min_probability * sum)) {
          *pmass *= (1. - nb_correction * min_probability) / norm;
        }
        reestim++;
        pmass++;
      }
    }
  }

  else {
    nb_correction = 0;
    for (i = 0;i < nb_value;i++) {
      if (*pmass++ > 0.) {
        nb_correction++;
      }
    }
    pmass -= nb_value;

    for (i = 0;i < nb_value;i++) {
      if (*pmass > 0.) {
        *pmass = 1. / (double)nb_correction;
      }
      pmass++;
    }
  }
}
