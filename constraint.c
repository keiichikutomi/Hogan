struct oconstraint{long int code,loff;
				   struct orevolute *revolutes;
				   struct ospherical *sphericals;
				   struct oprismatic *prismatics;
				   struct ocylindrical *cylindrical;


				  };



struct orevolute{long int code,loff;
				 struct onode *(node[2]);
				 double raxis[3];
				};/*1 DIRECTION ROTATION, NO TRANSITION.*/

struct ospherical{long int code,loff;
				  struct onode *(node[2]);
				 };/*FREE ROTATION, NO TRANSITION.*/

struct oprismatic{long int code,loff;
				  struct onode *(node[2]);
				  double saxis[3];
				 };/*1 DIRECTION TRANSITION, NO ROTATION.*/

struct ocylindrical{long int code,loff;
				   struct onode *(node[2]);
				   double saxis[3];
				   double raxis[3];
				  };/*1 DIRECTION ROTATION, 1 DIRECTION TRANSITION.*/

struct oplane{long int code,loff;
				  struct onode *(node[2]);

		     };

