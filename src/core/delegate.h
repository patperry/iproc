#ifndef _DELEGATE_H
#define _DELEGATE_H

typedef bool (*action_fn) (void *val, void *udata);
typedef bool (*predicate_fn) (const void *val, void *udata);

#endif /* _DELEGATE_H */
