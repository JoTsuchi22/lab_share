#ifndef PICK_H
#define PICK_H

/*
 * Handle node pick event
 */
//void handleNodePickEvent(void);

/*
 * Handle element pick event
 */
//void handleElementPickEvent(void);

void calculateStainStressIndexFromCoord(int patch_index, int index_xi, int index_eta,
                                        int *element_index_xi, int *element_index_eta,
                                        int *index_in_element_xi, int *index_in_element_eta,
                                        int *num_point_xi, int *num_point_eta);

/*
 * Handle point pick event
 */
void handlePointPickEvent(int patch_n, int *picked_point_flag);

#endif  /* PICK_H */
