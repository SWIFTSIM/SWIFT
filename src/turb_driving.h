

void do_turb_driving_step_first_half(struct cell *c, const struct engine *e);
void do_turb_driving_step_second_half(struct cell *c, const struct engine *e);

void init_turb(const struct space *s, const struct engine *e);
void add_turb_accel(const struct engine *e);

void log_turb_temp(struct engine *e);
